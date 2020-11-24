#include <igl/boundary_facets.h>
#include <igl/parallel_for.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/volume.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <prism/energy/smoother_pillar.hpp>
#include <prism/predicates/triangle_triangle_intersection.hpp>

#include "mesh_coloring.hpp"
#include "prism/PrismCage.hpp"
#include "prism/cage_utils.hpp"
#include "prism/cgal/triangle_triangle_intersection.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/intersections.hpp"
#include "prism/spatial-hash/AABB_hash.hpp"
#include "remesh_pass.hpp"
#include "validity_checks.hpp"

namespace prism::local_validity {
int attempt_relocate(PrismCage &pc, const std::vector<std::set<int>> &map_track,
                     const prism::local::RemeshOptions &option,
                     // specified infos below
                     const std::vector<int> &nb, int vid,
                     const std::array<Vec3d, 3> /*b-m-t*/ &relocations,
                     bool feature_enable,  // enable feature surgery for trackee
                     std::vector<std::set<int>> &sub_trackee) {
  auto &base = pc.base, &top = pc.top, &mid = pc.mid;
  auto &F = pc.F;
  auto &refV = pc.ref.V;
  auto &refF = pc.ref.F;
  auto &tree = *pc.ref.aabb;

  std::vector<Vec3i> nb_tris(nb.size());
  for (int i = 0; i < nb.size(); i++) nb_tris[i] = F[nb[i]];
  auto inter_check = [&pc, &option, &old_fid = nb, &moved_tris = nb_tris]() {
    if (!option.dynamic_hashgrid)
      return intersect_check(pc.base, pc.top, moved_tris, *pc.ref.aabb);
    else {
      const std::lock_guard<std::mutex> guard(pc.grid_mutex);
      return dynamic_intersect_check(pc.base, pc.F, old_fid, moved_tris,
                                     *pc.base_grid) &&
             dynamic_intersect_check(pc.top, pc.F, old_fid, moved_tris,
                                     *pc.top_grid);
    }
  };

  std::tuple<Vec3d, Vec3d, Vec3d> stored_loc = {base[vid], mid[vid], top[vid]};
  auto quality_before = max_quality_on_tris(base, mid, top, nb_tris);
  base[vid] = relocations[0];
  mid[vid] = relocations[1];
  top[vid] = relocations[2];
  auto after_failure = [&]() {
    std::tie(base[vid], mid[vid], top[vid]) = stored_loc;
  };

  auto quality_after = max_quality_on_tris(base, mid, top, nb_tris);
  if (quality_after > quality_before || std::isnan(quality_after)) {
    after_failure();
    return 4;
  }

  if (!volume_check(base, mid, top, nb_tris, tree.num_freeze)) {
    after_failure();
    return 1;
  }
  if (!inter_check()) {
    after_failure();
    spdlog::trace("failed inter check {}, restore", vid);
    return 2;
  }

  if (!feature_aware_distort_check(pc, map_track, option, nb, nb_tris, vid,
                                   feature_enable, sub_trackee)) {
    after_failure();
    spdlog::trace("failed map check {}, restore", vid);
    return 3;
  }
  spdlog::trace("Pass map check, moving on");
  return 0;
}
}  // namespace prism::local_validity

namespace prism::local {
void smooth_prism(PrismCage &pc, int vid,
                  const std::vector<std::vector<int>> &VF,
                  const std::vector<std::vector<int>> &VFi,
                  const RemeshOptions &option, bool true_zoom_false_rotate,
                  const std::vector<bool> &skip) {
  if (vid < pc.ref.aabb->num_freeze)
    return;  // only skip singularity, not boundary or feature
  std::optional<std::pair<Vec3d, Vec3d>> great_prism;

  if (true_zoom_false_rotate) {
    great_prism = zoom(pc.base, pc.mid, pc.top, pc.F, VF, VFi, vid,
                       option.target_thickness);
  } else
    great_prism = rotate(pc.base, pc.mid, pc.top, pc.F, VF, VFi, vid,
                         option.target_thickness);
  // great_prism =
  // zoom_and_rotate(pc.base, pc.mid, pc.top, pc.F, pc.ref.aabb->num_freeze, VF,
  // VFi, vid, target_thick);
  if (!great_prism) {
    spdlog::debug("Zoom and Rotate failed.");
    return;
  }

  std::array<Vec3d, 3> relocations;
  std::tie(relocations[0], relocations[2]) = great_prism.value();

  std::vector<std::set<int>> checker;
  auto flag = 1;
  double alpha = 1.;
  std::tuple<Vec3d, Vec3d, Vec3d> old_locations{pc.base[vid], pc.mid[vid],
                                                pc.top[vid]};
  bool enable_feature_separation =
      skip[vid];       // TODO, maybe split feature and real skip.
  while (flag == 1) {  // shrink if there is a volume failure.
    Vec3d new_b = relocations[0] * (alpha) + (1 - alpha) * pc.mid[vid];
    Vec3d new_t = relocations[2] * (alpha) + (1 - alpha) * pc.mid[vid];
    flag = prism::local_validity::attempt_relocate(
        pc, pc.track_ref, option, VF[vid], vid, {new_b, pc.mid[vid], new_t},
        enable_feature_separation, checker);
    spdlog::trace("newb {}", new_b);
    spdlog::trace("newt {}", new_t);
    alpha *= 0.8;
    if (alpha < 1e-2) break;
  }
  std::vector<RowMatd> local_cp;
  std::vector<Vec3i> moved_tris;
  if (flag == 0) {
    for (auto f : VF[vid]) moved_tris.push_back(pc.F[f]);
    if (option.curve_checker.first.has_value() &&
        !(std::any_cast<std::function<bool(
              const PrismCage &, const std::vector<int> &,
              const decltype(moved_tris) &, decltype(local_cp) &)>>(
            option.curve_checker.first))(pc, VF[vid], moved_tris, local_cp)) {
      flag = 5;
    }
  }
  if (flag > 0) {
    spdlog::debug("ZoomRotate checker failed.");
    std::tie(pc.base[vid], pc.mid[vid], pc.top[vid]) = old_locations;
    return;
  }
  if (option.curve_checker.second.has_value())
    (std::any_cast<
        std::function<void(const std::vector<int> &, const std::vector<int> &,
                           const std::vector<RowMatd> &)>>(
        option.curve_checker.second))(VF[vid], VF[vid], local_cp);

  if (pc.top_grid != nullptr) {
    spdlog::trace("HashGrid remove");
    const std::lock_guard<std::mutex> guard(pc.grid_mutex);
    for (auto f : VF[vid]) {
      pc.top_grid->remove_element(f);
      pc.base_grid->remove_element(f);
    }
    pc.top_grid->insert_triangles(pc.top, pc.F, VF[vid]);
    pc.base_grid->insert_triangles(pc.base, pc.F, VF[vid]);
  }
  for (int i = 0; i < VF[vid].size(); i++) {
    pc.track_ref[VF[vid][i]] = std::move(checker[i]);
  }
  spdlog::trace("ZoomRotate SUCCEED, move to next.");
}

void smooth_single(PrismCage &pc, int vid,
                   const std::vector<std::vector<int>> &VF,
                   const std::vector<std::vector<int>> &VFi,
                   const RemeshOptions &option, const std::vector<bool> &skip) {
  if (skip[vid]) return;

  spdlog::trace("smooth attempt: {}", vid);
  auto new_direction = prism::smoother_direction(
      pc.base, pc.mid, pc.top, pc.F, pc.ref.aabb->num_freeze, VF, VFi, vid);

  if (!new_direction) {
    spdlog::trace("No better location.");
    return;
  }

  std::tuple<Vec3d, Vec3d, Vec3d> old_locations{pc.base[vid], pc.mid[vid],
                                                pc.top[vid]};
  std::array<Vec3d, 3> relocations{pc.base[vid] + new_direction.value(),
                                   pc.mid[vid] + new_direction.value(),
                                   pc.top[vid] + new_direction.value()};
  auto query = [&ref = pc.ref](
                   const Vec3d &s, const Vec3d &t,
                   const std::set<int> &total_trackee) -> std::optional<Vec3d> {
    if (ref.aabb->enabled)  // this can be discarded if no performance benefit
                            // is found.
      return ref.aabb->segment_query(s, t);
    std::array<Vec3d, 2> seg_query{s, t};
    for (auto f : total_trackee) {
      auto v0 = ref.F(f, 0), v1 = ref.F(f, 1), v2 = ref.F(f, 2);
      auto mid_intersect = prism::cgal::segment_triangle_intersection(
          seg_query, {ref.V.row(v0), ref.V.row(v1), ref.V.row(v2)});
      if (mid_intersect) return mid_intersect;
    }
    return {};
  };

  if (true) {  // project onto reference
    std::set<int> total_trackee;
    for (auto f : VF[vid])
      total_trackee.insert(pc.track_ref[f].begin(), pc.track_ref[f].end());
    std::optional<Vec3d> mid_intersect;
    for (int i = 0; i < 20; i++) {
      mid_intersect = query(relocations[0], relocations[2], total_trackee);
      if (mid_intersect) break;
      relocations[0] = (pc.base[vid] + relocations[0]) / 2;
      relocations[2] = (pc.top[vid] + relocations[2]) / 2;
    }
    if (!mid_intersect) {
      spdlog::trace("Pan mid failed.");
      return;
    }
    relocations[1] = mid_intersect.value();
    spdlog::trace("found new mid {}", relocations[1]);
  }

  std::vector<std::set<int>> checker;
  int flag = prism::local_validity::attempt_relocate(
      pc, pc.track_ref, option, VF[vid], vid, relocations, {}, checker);

  std::vector<RowMatd> local_cp;
  std::vector<Vec3i> moved_tris;
  if (flag == 0) {
    for (auto f : VF[vid]) moved_tris.push_back(pc.F[f]);
    if (option.curve_checker.first.has_value() &&
        !(std::any_cast<std::function<bool(
              const PrismCage &, const std::vector<int> &,
              const decltype(moved_tris) &, decltype(local_cp) &)>>(
            option.curve_checker.first))(pc, VF[vid], moved_tris, local_cp)) {
      flag = 5;
    }
  }
  if (flag > 0) {
    spdlog::trace("Pan checker failed.");
    std::tie(pc.base[vid], pc.mid[vid], pc.top[vid]) = old_locations;
    return;
  }

 if (pc.top_grid != nullptr) {
    spdlog::trace("HashGrid remove");
    const std::lock_guard<std::mutex> guard(pc.grid_mutex);
    for (auto f : VF[vid]) {
      pc.top_grid->remove_element(f);
      pc.base_grid->remove_element(f);
    }
    pc.top_grid->insert_triangles(pc.top, pc.F, VF[vid]);
    pc.base_grid->insert_triangles(pc.base, pc.F, VF[vid]);
  }
  if (option.curve_checker.second.has_value())
    (std::any_cast<
        std::function<void(const std::vector<int> &, const std::vector<int> &,
                           const std::vector<RowMatd> &)>>(
        option.curve_checker.second))(VF[vid], VF[vid], local_cp);
  auto &new_tracks = checker;
  for (int i = 0; i < VF[vid].size(); i++) {
    pc.track_ref[VF[vid][i]] = std::move(new_tracks[i]);
  }
  spdlog::trace("SUCCEED, move to next.");
};
}  // namespace prism::local
void prism::local::localsmooth_pass(PrismCage &pc,
                                    const RemeshOptions &option) {
#ifndef NDEBUG
  {
    std::vector<Vec3d> tetV;
    std::vector<Vec4i> tetT;
    prism::cage_utils::tetmesh_from_prismcage(
        pc.base, pc.mid, pc.top, pc.F, pc.ref.aabb->num_freeze, tetV, tetT);
    RowMatd tV;
    RowMati tT;
    vec2eigen(tetV, tV);
    vec2eigen(tetT, tT);
    Eigen::VectorXd vols;
    igl::volume(tV, tT, vols);
    spdlog::warn("Volumes {} {} {}", vols.minCoeff(), vols.maxCoeff(),
                 vols.mean());
  }
#endif
  std::vector<std::vector<int>> VF, VFi, groups;
  std::vector<bool> skip_flag(pc.mid.size(), false);
  for (auto [m, ignore] : pc.meta_edges) {
    auto [v0, v1] = m;
    skip_flag[v0] = true;
    skip_flag[v1] = true;
  }
  for (int i = 0; i < pc.ref.aabb->num_freeze; i++) skip_flag[i] = true;
  {
    RowMati mF, mE;
    vec2eigen(pc.F, mF);
    igl::vertex_triangle_adjacency(pc.mid.size(), mF, VF, VFi);
    prism::local::vertex_coloring(mF, groups);
    igl::boundary_facets(mF, mE);
    spdlog::info("boundary, mE {}", mE.rows());
    for (int i = 0; i < mE.rows(); i++)
      for (auto j : {0, 1}) skip_flag[mE(i, j)] = true;
  }
  spdlog::info("Group Count {}", groups.size());

  std::srand(0);
  if (option.parallel && pc.top_grid != nullptr) {
    spdlog::error("Multithread hashmap is not safe. Todo: move to TBB.");
  }
  for (auto &gr : groups)
    igl::parallel_for(
        gr.size(),
        [&gr, &pc, &VF, &VFi, &option, &skip_flag](auto ii) {
          smooth_single(pc, gr[ii], VF, VFi, option, skip_flag);
        },
        size_t(option.parallel ? 1 : pc.mid.size()));
  for (auto &gr : groups)
    igl::parallel_for(
        gr.size(),
        [&gr, &pc, &VF, &VFi, &option, &skip_flag](auto ii) {
          smooth_prism(pc, gr[ii], VF, VFi, option, false, skip_flag);
          smooth_prism(pc, gr[ii], VF, VFi, option, true, skip_flag);
        },
        size_t(option.parallel ? 1 : pc.mid.size()));

  spdlog::info("Finished Smoothing");
  return;
}

namespace prism::local {

void legacy_smooth_prism(PrismCage &pc, int vid,
                         const std::vector<std::vector<int>> &VF,
                         const std::vector<std::vector<int>> &VFi,
                         const prism::local::RemeshOptions &option,
                         const std::vector<bool> &skip, bool on_base) {
  if (skip[vid]) return;

  auto great_prism = prism::smoother_location_legacy(
      pc.base, pc.mid, pc.top, pc.F, pc.ref.aabb->num_freeze, VF, VFi, vid,
      on_base);
  if (!great_prism) {
    spdlog::trace("Legacy Smooth failed.");
    return;
  }

  std::array<Vec3d, 3> relocations{pc.base[vid], pc.mid[vid], pc.top[vid]};
  relocations[on_base ? 0 : 2] = great_prism.value();
  relocations[1] = (relocations[0] + relocations[2]) / 2;

  std::vector<std::set<int>> checker;
  auto flag = 1;
  double alpha = 1.;
  flag = prism::local_validity::attempt_relocate(
      pc, pc.track_ref, option, VF[vid], vid, relocations, {}, checker);
  if (flag > 0) {
    spdlog::trace("ZoomRotate checker failed.");
  } else {
    for (int i = 0; i < VF[vid].size(); i++) {
      pc.track_ref[VF[vid][i]] = std::move(checker[i]);
    }
    spdlog::trace("ZoomRotate SUCCEED, move to next.");
  }
}

void shellsmooth_pass(PrismCage &pc, const RemeshOptions &option) {
  std::vector<std::vector<int>> VF, VFi, groups;
  std::vector<bool> skip_flag(pc.mid.size(), false);
  {
    RowMati mF, mE;
    vec2eigen(pc.F, mF);
    igl::vertex_triangle_adjacency(pc.mid.size(), mF, VF, VFi);
    prism::local::vertex_coloring(mF, groups);
    igl::boundary_facets(mF, mE);
    spdlog::info("boundary, mE {}", mE.rows());
    for (int i = 0; i < mE.rows(); i++)
      for (auto j : {0, 1}) {
        skip_flag[mE(i, j)] = true;
      }
    for (int i = 0; i < pc.ref.aabb->num_freeze; i++) skip_flag[i] = true;
  }
  spdlog::info("Group Count {}", groups.size());
  spdlog::set_level(spdlog::level::info);
  std::srand(0);

  for (auto &gr : groups)
    igl::parallel_for(
        gr.size(),
        [&gr, &pc, &VF, &VFi, &option, &skip_flag](auto ii) {
          legacy_smooth_prism(pc, gr[ii], VF, VFi, option, skip_flag, true);
          legacy_smooth_prism(pc, gr[ii], VF, VFi, option, skip_flag, false);
        },
        size_t(option.parallel ? 1 : pc.mid.size()));

  spdlog::set_level(spdlog::level::info);
  return;
}

}  // namespace prism::local
