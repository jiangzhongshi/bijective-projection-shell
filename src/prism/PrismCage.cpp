#include "PrismCage.hpp"

#include <igl/adjacency_list.h>
#include <igl/boundary_loop.h>
#include <igl/cat.h>
#include <igl/doublearea.h>
#include <igl/remove_unreferenced.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/volume.h>
#include <igl/writeOBJ.h>
#include <igl/write_triangle_mesh.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <highfive/H5Easy.hpp>
#include <prism/geogram/geogram_utils.hpp>
#include <prism/local_operations/validity_checks.hpp>
#include <stdexcept>

#include "bevel_utils.hpp"
#include "cage_utils.hpp"
#include "geogram/AABB.hpp"
#include "spatial-hash/AABB_hash.hpp"

PrismCage::PrismCage(const RowMatd &vert, const RowMati &face,
                     double doosabineps, double initial_step, SeparateType st) {
  prism::geo::init_geogram();
  RowMatd VN;
  ref.F = face;
  ref.V = vert;
  std::set<int> omni_singu;

  bool good_normals =
      prism::cage_utils::most_visible_normals(ref.V, ref.F, VN, omni_singu);
  int pure_singularity = omni_singu.size();
  prism::cage_utils::mark_singular_on_border(ref.V, ref.F, VN, omni_singu);

  if (omni_singu.empty()) {
    assert(good_normals);
    spdlog::info("Succeessful Normals");
  } else {
    spdlog::info("Omni Singularity {} (Border {})", omni_singu.size(),
                 omni_singu.size() - pure_singularity);
    spdlog::trace("<Freeze> Omni \n {}", (omni_singu));
  }

  int num_cons = omni_singu.size();

  std::set<int> feature_verts;
  // reorder vert to make omni in the front
  prism::cage_utils::reorder_singularity_to_front(
      ref.V, ref.F, VN, omni_singu, feature_verts, vertex_reorder);

  for (int i = 0; i < ref.F.rows(); i++) {
    auto [s, mt, shift] =
        tetra_split_AorB({ref.F(i, 0), ref.F(i, 1), ref.F(i, 2)});
    for (int j = 0; j < 3; j++) ref.F(i, j) = mt[j];
  }

  RowMatd dsV = ref.V, dsVN = VN;
  RowMati dsF = ref.F;
  Eigen::VectorXd M;
  igl::doublearea(ref.V, ref.F, M);
  spdlog::info("Area {}", M.minCoeff() / 2);
  std::vector<int> face_parent;
  prism::bevel_utils::adaptive_doo_sabin(ref.V, ref.F, VN, doosabineps, dsV,
                                         dsF, dsVN, face_parent);

  spdlog::info("V/F {}/{}-bevel-> {}/{} ", ref.V.rows(), ref.F.rows(),
               dsV.rows(), dsF.rows());

  if (dsVN.hasNaN()) exit(1);
  // special bevel for singularity.
  prism::bevel_utils::singularity_special_bevel(dsV, dsF, num_cons, dsVN,
                                                face_parent);

  for (int i = 0; i < dsF.rows(); i++) {
    auto [s, mt, shift] = tetra_split_AorB({dsF(i, 0), dsF(i, 1), dsF(i, 2)});
    for (int j = 0; j < 3; j++) dsF(i, j) = mt[j];
  }

  ref.aabb.reset(
      new prism::geogram::AABB(ref.V, ref.F, st == SeparateType::kSurface));
  ref.aabb->num_freeze = num_cons;

  RowMatd inner, outer;
  prism::cage_utils::extrude_for_base_and_top(
      dsV, dsF, *ref.aabb, dsVN, num_cons, inner, outer, initial_step);
  eigen2vec(outer, top);
  eigen2vec(inner, base);

  eigen2vec(dsV, mid);
  eigen2vec(dsF, F);

  std::vector<std::vector<int>> VF, VFi;
  igl::vertex_triangle_adjacency(dsV, dsF, VF, VFi);
  if (st == SeparateType::kShell) {
    prism::cage_utils::hashgrid_shrink(mid, top, F, VF);
    prism::cage_utils::hashgrid_shrink(mid, base, F, VF);
    spdlog::warn(
        "Using hashgrid collision detection. Improvement in progress.");
    top_grid.reset(new prism::HashGrid(top, F));
    base_grid.reset(new prism::HashGrid(base, F));
  }
  // initial track with itself
  track_ref.resize(F.size());

  for (auto i = 0; i < F.size(); i++) {
    for (auto j : {0, 1, 2}) {
      auto v = F[i][j];
      if (v < num_cons) continue;
      for (auto t : VF[v]) {
        track_ref[i].insert(face_parent[t]);
      }
    }
  }
}

void PrismCage::init_track() {
  if (mid.size() == 0) eigen2vec(ref.V, mid);
  if (F.size() == 0) eigen2vec(ref.F, F);

#ifndef NDEBUG
  // Reorder F
  for (int i = 0; i < F.size(); i++) {
    auto [s, mt, shift] = tetra_split_AorB(F[i]);
    F[i] = mt;
    assert(shift == 0);
  }
#endif
  // initial track with itself
  assert(track_ref.size() == F.size());
}

auto serialize_meta_edges = [](auto &meta_edges) {
  std::vector<int> flat, ind;
  for (auto [m, data] : meta_edges) {
    ind.push_back(flat.size());
    flat.push_back(m.first);
    flat.push_back(m.second);
    auto [cid, seg] = data;
    flat.push_back(cid);
    flat.insert(flat.end(), seg.begin(), seg.end());
  }
  ind.push_back(flat.size());
  return std::tuple(flat, ind);
};

auto deserialize_meta_edges = [](auto &flat, auto &ind) {
  std::map<std::pair<int, int>, std::pair<int, std::vector<int>>> meta;
  auto prev = 0;
  auto begin = flat.begin();
  for (auto cur : ind) {
    if (cur == 0) continue;  // compatible with 0 leading or not.
    // [prev, cur)
    auto itv0 = begin + prev;
    auto itv1 = std::next(itv0);
    std::vector<int> vec;
    std::copy(itv0 + 2 + 1, begin + cur, std::back_inserter(vec));
    meta[{*itv0, *itv1}] = {*(itv0 + 2), vec};
    prev = cur;
  }
  return meta;
};

void PrismCage::serialize(std::string filename, std::any additionals) {
  RowMatd mbase, mtop, mV;
  RowMati mF;
  vec2eigen(base, mbase);
  vec2eigen(top, mtop);
  vec2eigen(mid, mV);
  vec2eigen(F, mF);

  Eigen::RowVectorXi track_sizes(track_ref.size());
  std::vector<int> track_flat;
  for (int i = 0; i < track_ref.size(); i++) {
    track_sizes[i] = track_ref[i].size();
    for (auto t : track_ref[i]) {
      track_flat.push_back(t);
    }
  }

  H5Easy::File file(filename, H5Easy::File::Overwrite);

  H5Easy::dump(file, "ref.V", ref.V);
  H5Easy::dump(file, "ref.F", ref.F);
  H5Easy::dump(file, "mbase", mbase);
  H5Easy::dump(file, "mtop", mtop);
  H5Easy::dump(file, "mV", mV);
  H5Easy::dump(file, "mF", mF);
  H5Easy::dump(file, "track_flat", track_flat);
  H5Easy::dump(file, "track_size", track_sizes);
  auto [meta_edges_flat, meta_edges_ind] = serialize_meta_edges(meta_edges);
  H5Easy::dump(file, "meta_edges_flat", meta_edges_flat);
  H5Easy::dump(file, "meta_edges_ind", meta_edges_ind);

  SeparateType st =
      ref.aabb->enabled ? SeparateType::kSurface : SeparateType::kNone;
  if (top_grid != nullptr) st = SeparateType::kShell;
  H5Easy::dump(file, "metadata",
               std::vector<int>{1, static_cast<int>(st), ref.aabb->num_freeze});
  spdlog::info("SeparateType {}",
               st == SeparateType::kSurface
                   ? "kSurface"
                   : (st == SeparateType::kShell ? "kShell" : "kNone"));

  if (additionals.has_value())
    std::any_cast<std::function<void(decltype(file) &)>>(additionals)(file);
}

void PrismCage::load_from_hdf5(std::string filename) {
  H5Easy::File file(filename, H5Easy::File::ReadOnly);
  RowMatd mbase, mtop, mV;
  RowMati mF;
  Eigen::RowVectorXi mTracks;
  Eigen::RowVectorXi track_size;
  std::vector<int> metadata;
  ref.V = H5Easy::load<decltype(ref.V)>(file, "ref.V");
  ref.F = H5Easy::load<decltype(ref.F)>(file, "ref.F");
  mbase = H5Easy::load<decltype(mbase)>(file, "mbase");
  mtop = H5Easy::load<decltype(mtop)>(file, "mtop");
  mV = H5Easy::load<decltype(mV)>(file, "mV");
  mF = H5Easy::load<decltype(mF)>(file, "mF");
  mTracks = H5Easy::load<decltype(mTracks)>(file, "track_flat");
  track_size = H5Easy::load<decltype(track_size)>(file, "track_size");
  SeparateType st = SeparateType::kSurface;
  if (file.exist("metadata")) {
    metadata = H5Easy::load<decltype(metadata)>(file, "metadata");
    st = static_cast<SeparateType>(metadata[1]);
  }
  if (file.exist("meta_edges_flat")) {
    std::vector<int> flat, ind;
    flat = H5Easy::load<decltype(flat)>(file, "meta_edges_flat");
    ind = H5Easy::load<decltype(ind)>(file, "meta_edges_ind");
    spdlog::debug("loading meta edge, with info {} {}", flat.size(),
                  ind.size());
    meta_edges = deserialize_meta_edges(flat, ind);
  }

  eigen2vec(mbase, base);
  eigen2vec(mtop, top);
  eigen2vec(mV, mid);
  eigen2vec(mF, F);

  // after
  ref.aabb = std::make_unique<prism::geogram::AABB>(
      ref.V, ref.F, st == SeparateType::kSurface);
  for (int i = 0; i < mid.size(); i++) {
    if (base[i] != top[i]) {
      ref.aabb->num_freeze = i;
      break;
    }
  }
  spdlog::info("AABB {}, SeparateType {}",
               ref.aabb->enabled ? "enabled" : "disabled",
               st == SeparateType::kSurface
                   ? "kSurface"
                   : (st == SeparateType::kShell ? "kShell" : "kNone"));
  if (st == SeparateType::kShell) {
    spdlog::info("Loading HashGrid.");
    top_grid.reset(new prism::HashGrid(top, F));
    base_grid.reset(new prism::HashGrid(base, F));
  }
  spdlog::info("Loaded singularity {}", ref.aabb->num_freeze);
  track_ref.clear();
  for (int i = 0, cur = 0; i < track_size.size(); i++) {
    auto ts = track_size[i];
    std::set<int> cur_track(mTracks.data() + cur, mTracks.data() + cur + ts);
    track_ref.emplace_back(cur_track);
    cur += ts;
  }
}

PrismCage::PrismCage(std::string filename) {
  prism::geo::init_geogram();
  namespace fs = std::filesystem;
  assert(fs::path(filename).extension() == ".h5" ||
         fs::path(filename).extension() == ".init");
  spdlog::info("Loading From .H5");
  load_from_hdf5(filename);
}

void PrismCage::cleanup_empty_faces(Eigen::VectorXi &NI, Eigen::VectorXi &NJ) {
  // this is called after collapse pass.
  constexpr auto mark_zero_rows = [](const auto &vecF, RowMati &mat) {
    std::vector<Vec3i> newF;
    newF.reserve(vecF.size());
    for (int i = 0; i < vecF.size(); i++) {
      auto &f = vecF[i];
      if (f[0] != f[1])
        newF.push_back(f);
      else
        newF.push_back({-1, -1, -1});
    }
    mat = Eigen::Map<RowMati>(newF[0].data(), newF.size(), 3);
  };

  RowMati mF;
  mark_zero_rows(F, mF);
  igl::remove_unreferenced(mid.size(), mF, NI, NJ);

#ifndef NDEBUG
  // assuming NJ is sorted ascending
  for (int i = 0; i < NJ.size() - 1; i++) assert(NJ[i] < NJ[i + 1]);
#endif
  for (int i = 0; i < NJ.size(); i++) {
    mid[i] = mid[NJ[i]];
    base[i] = base[NJ[i]];
    top[i] = top[NJ[i]];
  }
  mid.resize(NJ.size());
  base.resize(NJ.size());
  top.resize(NJ.size());

  int cur = 0;
  for (int i = 0; i < F.size(); i++) {
    if (F[i][0] == F[i][1]) continue;
    if (track_ref[i].size() == 0) spdlog::error("Zero Tracer");
    if (i != cur) track_ref[cur] = std::move(track_ref[i]);
    for (int j = 0; j < 3; j++) F[cur][j] = NI[F[i][j]];
    if (F[cur][0] > F[cur][1] || F[cur][0] > F[cur][2])
      spdlog::error("v0 v1 v2 order wrong at {}", cur);
    cur++;
  }
  track_ref.resize(cur);
  F.resize(cur);

  auto &vid_map = NI;
  // feature meta edges
  std::map<std::pair<int, int>, std::pair<int, std::vector<int>>> new_metas;
  for (auto m : meta_edges) {
    auto [u0, u1] = m.first;
    new_metas[{vid_map[u0], vid_map[u1]}] = m.second;
  }
  meta_edges = std::move(new_metas);

  // hash grid update.
  if (top_grid != nullptr) {
    top_grid->update_after_collapse();
    base_grid->update_after_collapse();
    assert(top_grid->face_stores.size() == F.size());
    assert(base_grid->face_stores.size() == F.size());
  }
}