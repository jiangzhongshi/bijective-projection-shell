#include "validity_checks.hpp"

#include <igl/Timer.h>
#include <igl/parallel_for.h>
#include <igl/volume.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <highfive/H5Easy.hpp>
#include <limits>
#include <prism/predicates/inside_octahedron.hpp>
#include <prism/predicates/triangle_triangle_intersection.hpp>
#include <queue>
#include <vector>

#include "../energy/map_distortion.hpp"
#include "../energy/prism_quality.hpp"
#include "../geogram/AABB.hpp"
#include "../predicates/inside_prism_tetra.hpp"
#include "../predicates/positive_prism_volume_12.hpp"
#include "../spatial-hash/AABB_hash.hpp"
#include "local_mesh_edit.hpp"
#include "prism/PrismCage.hpp"
#include "remesh_pass.hpp"

namespace prism::local_validity {

constexpr auto quality_on_tris = [](const auto &base, const auto &mid,
                                    const auto &top, const auto &moved_tris) {
  double quality = 0;

  for (auto [v0, v1, v2] : moved_tris) {
    auto q = prism::energy::triangle_quality({mid[v0], mid[v1], mid[v2]});
    quality += q;
  }
  return quality;
};

double max_quality_on_tris(const std::vector<Vec3d> &base,
                           const std::vector<Vec3d> &mid,
                           const std::vector<Vec3d> &top,
                           const std::vector<Vec3i> &moved_tris) {
  double quality = 0;

  for (auto [v0, v1, v2] : moved_tris) {
    // auto q = prism::energy::prism_full_quality(
    //              {base[v0], base[v1], base[v2], mid[v0], mid[v1], mid[v2]}) +
    //          prism::energy::prism_full_quality(
    //              {mid[v0], mid[v1], mid[v2], top[v0], top[v1], top[v2]});
    auto q = prism::energy::triangle_quality({mid[v0], mid[v1], mid[v2]});
    if (std::isnan(q)) return std::numeric_limits<double>::infinity();
    quality = std::max(quality, q);
  }
  return quality;
}

bool dynamic_intersect_check(
    const std::vector<Vec3d> &base, const std::vector<Vec3i> &F,
    const std::vector<int>
        &vec_removed,  // proposed removal face_id to be ignored in the test.
    const std::vector<Vec3i> &tris,  // proposed addition triangles
    const prism::HashGrid &grid) {
  spdlog::trace("In DIC 2x{}", tris.size());
  std::set<int> removed(
      vec_removed.begin(),
      vec_removed.end());  // important, this has to be sorted or set for
                           // std::difference to work
  constexpr auto share_vertex = [](const auto &F, int c, auto &f) {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (F[c][i] == f[j]) return true;
    return false;  // TODO prove this complete, one vertex touch requires
                   // reduced collision check.
  };
  constexpr auto intersection = [](const auto &V, const auto &F, int c,
                                   auto &f) {
    auto &[p0, p1, p2] = F[c];
    auto &[q0, q1, q2] = f;
    return prism::predicates::triangle_triangle_overlap({V[p0], V[p1], V[p2]},
                                                        {V[q0], V[q1], V[q2]});
  };
  igl::Timer timer;
  timer.start();
  for (auto f : tris) {
    RowMat3d local;
    for (auto k : {0, 1, 2}) local.row(k) = base[f[k]];
    std::set<int> candidates;
    grid.query(local.colwise().minCoeff(), local.colwise().maxCoeff(),
               candidates);
    std::vector<int> result;
    std::set_difference(candidates.begin(), candidates.end(), removed.begin(),
                        removed.end(), std::back_inserter(result));
    if (result.empty()) {
      spdlog::error("Empty query for spatial hash.");
    }
    for (auto c : result) {  // candidate faces to test.
      if ((!share_vertex(F, c, f)) && intersection(base, F, c, f)) return false;
    }
  }
  auto elapsed = timer.getElapsedTimeInMicroSec();
  spdlog::trace("DIC true {}", elapsed);
  return true;  // safe operation, no intersection
}

bool prism_positivity_with_numerical(const std::array<Vec3d, 6> &verts,
                                     const std::array<bool, 3> &constrained) {
  auto matV = Eigen::Map<const RowMatd>(verts[0].data(), 6, 3);
  for (int p = 0; p < 3; p++) {
    if (constrained[p]) continue;
    Eigen::VectorXd vol_index;
    igl::volume(matV,
                Eigen::Map<const RowMati>(TWELVE_TETRAS[4 * p].data(), 4, 4),
                vol_index);
    if (vol_index.minCoeff() <= 0) {
      spdlog::trace(vol_index.format(Eigen::IOFormat(Eigen::FullPrecision)));
      return false;
    }
  }
  return prism::predicates::positive_prism_volume(verts, constrained);
};

bool volume_check(const std::vector<Vec3d> &base, const std::vector<Vec3d> &mid,
                  const std::vector<Vec3d> &top, const std::vector<Vec3i> &tris,
                  int num_cons) {
  //
  spdlog::trace("In VC");
  igl::Timer timer;
  timer.start();
  for (auto [v0, v1, v2] : tris) {
    std::array<bool, 3> cons_flag{v0 < num_cons, v1 < num_cons, v2 < num_cons};
    if (!prism::predicates::positive_prism_volume(
            {base[v0], base[v1], base[v2], mid[v0], mid[v1], mid[v2]},
            cons_flag))
      return false;
    if (!prism::predicates::positive_prism_volume(
            {mid[v0], mid[v1], mid[v2], top[v0], top[v1], top[v2]}, cons_flag))
      return false;
  }
  auto elapsed = timer.getElapsedTimeInMicroSec();
  spdlog::trace("VC true {}", elapsed);
  return true;
}

bool volume_check(const std::vector<Vec3d> &base, const std::vector<Vec3d> &top,
                  const std::vector<Vec3i> &tris, int num_cons) {
  //
  spdlog::trace("In VC");
  igl::Timer timer;
  timer.start();
  for (auto [v0, v1, v2] : tris) {
    std::array<bool, 3> cons_flag{v0 < num_cons, v1 < num_cons, v2 < num_cons};
    if (!prism::predicates::positive_prism_volume(
            {base[v0], base[v1], base[v2], top[v0], top[v1], top[v2]},
            cons_flag))
      return false;
  }
  auto elapsed = timer.getElapsedTimeInMicroSec();
  spdlog::trace("VC true {}", elapsed);
  return true;
}

bool intersect_check(const std::vector<Vec3d> &base,
                     const std::vector<Vec3d> &top,
                     const std::vector<Vec3i> &tris,
                     const prism::geogram::AABB &tree) {
  spdlog::trace("In IC 2x{}", tris.size());
  igl::Timer timer;
  timer.start();
  for (auto [v0, v1, v2] : tris) {
    if (tree.intersects_triangle({base[v0], base[v1], base[v2]},
                                 v0 < tree.num_freeze)) {
      spdlog::trace("base {} {} {}", v0, v1, v2);
      return false;
    }
    if (tree.intersects_triangle({top[v0], top[v1], top[v2]},
                                 v0 < tree.num_freeze)) {
      spdlog::trace("top {} {} {}", v0, v1, v2);
      return false;
    }
  }
  auto elapsed = timer.getElapsedTimeInMicroSec();
  spdlog::trace("IC true {}", elapsed);
  return true;
}

std::optional<std::vector<std::set<int>>> distort_check(
    const std::vector<Vec3d> &base,
    const std::vector<Vec3d> &mid,  // placed new verts
    const std::vector<Vec3d> &top, const std::vector<Vec3i> &tris,
    const std::set<int> &combined_trackee,  // indices to ref.F tracked
    const RowMatd &refV, const RowMati &refF, double distortion_bound,
    int num_freeze, bool bundled_intersection) {
  // NormalCheck
  spdlog::trace("In NC ct#{}, tris{}", combined_trackee.size(), tris.size());
  igl::Timer timer;
  timer.start();
  assert(base.size() == top.size());
  assert(base.size() == mid.size());
  std::vector<std::set<int>> distributed_refs(tris.size());
  for (int i = 0; i < tris.size(); i++) {
    auto [v0, v1, v2] = tris[i];
    std::array<Vec3d, 3> base_vert{base[v0], base[v1], base[v2]};
    std::array<Vec3d, 3> mid_vert{mid[v0], mid[v1], mid[v2]};
    std::array<Vec3d, 3> top_vert{top[v0], top[v1], top[v2]};
    std::array<bool, 3> oct_type_b, oct_type_t;
    prism::determine_convex_octahedron(base_vert, mid_vert, oct_type_b,
                                       num_freeze > v0);
    prism::determine_convex_octahedron(mid_vert, top_vert, oct_type_t,
                                       num_freeze > v0);
    spdlog::trace("checking tris{}: {}-{}-{}", i, v0, v1, v2);
    for (auto t : combined_trackee) {  // for every tracked original triangle.
      std::array<Vec3d, 3> ref_tri = {
          refV.row(refF(t, 0)), refV.row(refF(t, 1)), refV.row(refF(t, 2))};
      bool intersected_prism = false;
      if (v0 >= num_freeze || v0 != refF(t, 0)) {
        intersected_prism =
            prism::triangle_intersect_octahedron(
                base_vert, mid_vert, oct_type_b, ref_tri, num_freeze > v0) ||
            prism::triangle_intersect_octahedron(mid_vert, top_vert, oct_type_t,
                                                 ref_tri, num_freeze > v0);
      } else
        intersected_prism = prism::singularless_triangle_intersect_octahedron(
                                base_vert, mid_vert, oct_type_b, ref_tri) ||
                            prism::singularless_triangle_intersect_octahedron(
                                mid_vert, top_vert, oct_type_t, ref_tri);
      if (!intersected_prism) continue;

      for (int tc = (v0 < num_freeze) ? 1 : 0; tc < 3; tc++) {
        auto pillar = top_vert[tc] - base_vert[tc];
        auto distortion = prism::energy::map_max_cos_angle(pillar, ref_tri);
        if (distortion < distortion_bound) {
          spdlog::trace("ref{} tris {}-{}-{}, tc{}, distortion: {}", t, v0, v1,
                        v2, tc, distortion);
          return {};
        }
      }
      distributed_refs[i].insert(t);
      // Test intersect between tri(t) and top/base(i)
      if (bundled_intersection &&  // this is enabled when global (AABB tree)
                                   // intersection is off.
          (prism::predicates::triangle_triangle_overlap(base_vert, ref_tri) ||
           prism::predicates::triangle_triangle_overlap(top_vert, ref_tri))) {
        spdlog::trace("bundled overlap detected");
        return {};
      }
    }
  }
  spdlog::trace("distributed {}", distributed_refs);

  auto elapsed = timer.getElapsedTimeInMicroSec();
  spdlog::trace("NC true {}", elapsed);

  //
  return distributed_refs;
}

auto trackee_based_component_color = [](auto &map_track, auto &nb) {
  auto N = nb.size();
  std::vector<std::vector<int>> adj_list(N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) {
      if (non_empty_intersect(map_track[nb[i]], map_track[nb[j]])) {
        adj_list[i].push_back(j);
        adj_list[j].push_back(i);
      }
    }
  }
  std::vector<bool> visited(N, false);
  std::deque<int> bfsq;
  std::vector<int> separate_type(nb.size(), -1);
  auto current_color = 0;
  while (true) {
    auto nz_it = std::find(separate_type.begin(), separate_type.end(), -1);
    if (nz_it == separate_type.end()) break;
    bfsq.push_back(std::distance(separate_type.begin(), nz_it));
    do {
      auto id = bfsq.front();
      bfsq.pop_front();
      separate_type[id] = current_color;
      for (auto n : adj_list[id])
        if (separate_type[n] == -1) bfsq.push_back(n);
    } while (!bfsq.empty());
    current_color++;
  }
  return separate_type;
};

bool special_handling_hanging_terminal(
    const PrismCage &pc, const std::vector<std::set<int>> &map_track,
    const prism::local::RemeshOptions &option, const std::vector<int> &nb,
    const std::vector<Vec3i> &nb_tris, int vid,
    std::vector<std::set<int>> &sub_trackee) {
  auto &base = pc.base, &top = pc.top, &mid = pc.mid;
  auto &F = pc.F;
  auto &refV = pc.ref.V;
  auto &refF = pc.ref.F;
  auto &tree = *pc.ref.aabb;
  auto set_minus = [](const auto &A, const auto &B, auto &C) {
    std::set_difference(A.begin(), A.end(), B.begin(), B.end(),
                        std::inserter(C, C.begin()));
  };
  std::set<int> combined_tracks;
  for (auto f : nb) set_add_to(map_track[f], combined_tracks);

  std::vector<Vec3i> ordinary_checks;
  spdlog::trace("nbtris {}", nb_tris);
  for (auto f : nb_tris) {
    std::set<int> trackee;
    int j = 0;
    for (; j < 3; j++)
      if (f[j] == vid) break;
    assert(j != 3);
    auto v0 = vid, v1 = f[(j + 1) % 3], v2 = f[(j + 2) % 3];
    spdlog::trace("v0 {} v1 {} v2 {}", v0, v1, v2);
    auto it0 = pc.meta_edges.find({v0, v1});
    if (it0 == pc.meta_edges.end()) it0 = pc.meta_edges.find({v2, v0});
    auto it1 = pc.meta_edges.find({v1, v0});
    if (it1 == pc.meta_edges.end()) it1 = pc.meta_edges.find({v0, v2});
    if (it0 == it1) {  // normal check
      trackee = combined_tracks;
    } else {
      int swapped = 0;
      if (it0 == pc.meta_edges.end()) {
        swapped = 1;
        it0 = it1;
      }
      auto [cid, segs] = it0->second;
      set_minus(combined_tracks, option.chain_reject_trackee[cid * 2 + swapped],
                trackee);
      spdlog::trace("REJ {}", option.chain_reject_trackee[cid * 2 + swapped]);
      spdlog::trace("cid {}, swap {} trackee {}", cid, swapped, trackee);
    }
    auto sub_refs = distort_check(base, mid, top, {f}, trackee, refV, refF,
                                  option.distortion_bound, tree.num_freeze,
                                  option.dynamic_hashgrid);
    if (!sub_refs) return false;
    sub_trackee.emplace_back(sub_refs.value()[0]);
  }
  return true;
}

auto one_ring_component_coloring = [](const PrismCage &pc,
                                      const std::vector<int> &nb,
                                      const std::vector<Vec3i> &nb_tris,
                                      int vid) {
  std::map<int, std::array<int, 2>> vert2face;
  auto found_in = [](auto &map, auto &v0, auto &v1) {
    return map.find({v0, v1}) != map.end() || map.find({v1, v0}) != map.end();
  };
  for (int fi = 0; fi < nb_tris.size(); fi++) {
    auto &f = nb_tris[fi];
    std::set<int> trackee;
    int j = 0;
    for (; j < 3; j++)
      if (f[j] == vid) break;
    assert(j != 3);
    auto v0 = vid, v1 = f[(j + 1) % 3], v2 = f[(j + 2) % 3];

    for (auto vv : {v1, v2}) {
      if (found_in(pc.meta_edges, v0, vv)) {
        spdlog::trace("feat {} {}", v0, vv);
        continue;
      }
      auto it = vert2face.lower_bound(vv);
      if (it != vert2face.end() && it->first == vv)
        it->second[1] = (fi);
      else
        vert2face.emplace_hint(it, vv, std::array<int, 2>{fi, -1});
    }
  }
  spdlog::trace("nb_tris {}", nb_tris);
  spdlog::trace("v2f {}", vert2face);

  std::vector<int> separate_type(nb_tris.size(), -1);
  int color = 0;
  for (color = 0; !vert2face.empty(); color++) {
    std::deque<int> fid_queue;
    fid_queue.push_back(vert2face.begin()->first);
    while (!fid_queue.empty()) {
      spdlog::trace("fq {}color {}", fid_queue.size(), color);
      auto v1 = fid_queue.front();
      fid_queue.pop_front();
      auto it = vert2face.find(v1);
      if (it != vert2face.end()) {
        for (auto f : it->second) {
          separate_type[f] = color;
          for (auto j = 0; j < 3; j++) {
            auto v2 = nb_tris[f][j];
            if (v2 != v1 && v2 != vid) fid_queue.push_back(v2);
          }
        }
        vert2face.erase(it);
      }
    }
    if (color > 2 * nb_tris.size()) spdlog::error("inf loop");
  }
  for (auto &s : separate_type) {
    if (s == -1) s = color++;
  }
  spdlog::trace("sep {}", separate_type);
  return separate_type;
};

bool feature_aware_distort_check(const PrismCage &pc,
                                 const std::vector<std::set<int>> &map_track,
                                 const prism::local::RemeshOptions &option,
                                 const std::vector<int> &nb,
                                 const std::vector<Vec3i> &nb_tris, int vid,
                                 bool feature_enable,
                                 std::vector<std::set<int>> &sub_trackee) {
  auto &base = pc.base, &top = pc.top, &mid = pc.mid;
  auto &F = pc.F;
  auto &refV = pc.ref.V;
  auto &refF = pc.ref.F;
  auto &tree = *pc.ref.aabb;

  if (feature_enable) {  // special surgery.
    spdlog::trace("vid {}", vid);
    auto separate_type = one_ring_component_coloring(pc, nb, nb_tris, vid);
    spdlog::trace(separate_type);
    auto num_class =
        *std::max_element(separate_type.begin(), separate_type.end()) + 1;
    if (num_class == 1) {
      bool special_pass = special_handling_hanging_terminal(
          pc, map_track, option, nb, nb_tris, vid, sub_trackee);
      spdlog::trace("hanging terminal smoothing {}.", special_pass);
      return special_pass;
    }
    std::vector<std::vector<Vec3i>> moved_all(num_class);
    std::vector<std::set<int>> side_combined(num_class);
    for (int i = 0; i < separate_type.size(); i++) {
      auto f = nb[i];
      auto id = separate_type[i];
      moved_all[id].push_back(nb_tris[i]);
      set_add_to(map_track[f], side_combined[id]);
    }
    std::vector<std::vector<std::set<int>>> all_sub_refs;
    for (int i = 0; i < moved_all.size(); i++) {
      auto sub_refs_f = distort_check(
          base, mid, top, moved_all[i], side_combined[i], refV, refF,
          option.distortion_bound, tree.num_freeze, option.dynamic_hashgrid);
      if (!sub_refs_f) return false;
      all_sub_refs.emplace_back(sub_refs_f.value());
    }
    assert(separate_type.size() == nb.size());
    std::vector<int> l(num_class, 0);
    for (auto c : separate_type) {
      sub_trackee.emplace_back(std::move(all_sub_refs[c][l[c]++]));
    }
    spdlog::debug("Sucessfully perform a surgical relocate check");
  } else {  // ordinary relocate and check.
    std::set<int> combined_tracks;
    for (auto f : nb) {
      std::merge(map_track[f].begin(), map_track[f].end(),
                 combined_tracks.begin(), combined_tracks.end(),
                 std::inserter(combined_tracks, combined_tracks.begin()));
    }
    auto sub_refs = distort_check(base, mid, top, nb_tris, combined_tracks,
                                  refV, refF, option.distortion_bound,
                                  tree.num_freeze, option.dynamic_hashgrid);
    if (!sub_refs) {
      return false;
    }
    spdlog::trace("Pass map check, moving on");
    sub_trackee = std::move(sub_refs.value());
  }
  return true;
}

// below is the checks for local operations

int attempt_local_edit(
    const PrismCage &pc, const std::vector<std::set<int>> &map_track,
    const prism::local::RemeshOptions &option, const std::vector<int> &old_fid,
    const std::vector<Vec3i> &old_tris,
    const std::vector<Vec3i> &moved_tris,  // should be already shifted.
    std::tuple<std::vector<int> /*newfid*/, std::vector<int> /*shifts*/,
               std::vector<std::set<int>> /*track*/> &checker) {
  auto &base = pc.base, &top = pc.top, &mid = pc.mid;
  auto &F = pc.F;
  auto &refV = pc.ref.V;
  auto &refF = pc.ref.F;
  auto &tree = *pc.ref.aabb;

  std::vector<int> new_fid;

  auto inter_safe = [&pc, &option, &old_fid, &moved_tris]() {
    if (!option.dynamic_hashgrid)
      return intersect_check(pc.base, pc.top, moved_tris, *pc.ref.aabb);
    else
      return dynamic_intersect_check(pc.base, pc.F, old_fid, moved_tris,
                                     *pc.base_grid) &&
             dynamic_intersect_check(pc.top, pc.F, old_fid, moved_tris,
                                     *pc.top_grid);
  };
  // note that new_fid is geometric cyclic.

  spdlog::trace("Quality check");
  auto quality_before = max_quality_on_tris(base, mid, top, old_tris);
  auto quality_after = max_quality_on_tris(base, mid, top, moved_tris);
  spdlog::trace("Quality compare {} -> {}", quality_before, quality_after);
  if (std::isnan(quality_after)) return 4;
  if ((quality_after > option.collapse_quality_threshold) &&
      quality_after > quality_before)  // if the quality is too large, not allow
                                       // it to increase.
    return 4;
  //  volume check
  auto vc = volume_check(base, mid, top, moved_tris, tree.num_freeze);
  if (!vc) return 1;

  if (!inter_safe()) return 2;
  // get new subdivide types
  // auto new_shifts = prism::local_validity::triangle_shifts(moved_tris);

  std::set<int> combined_tracks;
  for (auto f : old_fid) {
    std::merge(map_track[f].begin(), map_track[f].end(),
               combined_tracks.begin(), combined_tracks.end(),
               std::inserter(combined_tracks, combined_tracks.begin()));
  }
  auto sub_refs = distort_check(base, mid, top, moved_tris, combined_tracks,
                                refV, refF, option.distortion_bound,
                                tree.num_freeze, option.dynamic_hashgrid);
  if (!sub_refs) return 3;

  // insert curve_check here.

  checker = std::tuple(std::move(new_fid), std::move(std::vector<int>()),
                       std::move(sub_refs.value()));
  return 0;
}

}  // namespace prism::local_validity