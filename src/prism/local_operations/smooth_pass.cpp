#include "mesh_coloring.hpp"
#include "prism/PrismCage.hpp"
#include "prism/cage_utils.hpp"
#include "prism/cgal/triangle_triangle_intersection.hpp"
#include "prism/geogram/AABB.hpp"
#include "remesh_pass.hpp"
#include "validity_checks.hpp"
#include <algorithm>
#include <igl/boundary_facets.h>
#include <igl/parallel_for.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/volume.h>
#include <prism/energy/smoother_pillar.hpp>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>
#include "prism/spatial-hash/AABB_hash.hpp"

namespace prism::local {
void smooth_prism(PrismCage &pc, int vid,
                  const std::vector<std::vector<int>> &VF,
                  const std::vector<std::vector<int>> &VFi,
                  double distortion_bound, double target_thick,
                  bool true_zoom_false_rotate, const std::vector<bool> &skip) {
  if (vid < pc.ref.aabb->num_freeze)
    return; // only skip singularity, not boundary or feature
  if (skip[vid])
    return; // TODO, deal with special distortion check for feature verts.
  std::optional<std::pair<Vec3d, Vec3d>> great_prism;

  if (true_zoom_false_rotate) {
    great_prism =
        zoom(pc.base, pc.mid, pc.top, pc.F, VF, VFi, vid, target_thick);
  } else
    great_prism =
        rotate(pc.base, pc.mid, pc.top, pc.F, VF, VFi, vid, target_thick);
  // great_prism =
  // zoom_and_rotate(pc.base, pc.mid, pc.top, pc.F, pc.ref.aabb->num_freeze, VF,
  // VFi, vid, target_thick);
  if (!great_prism) {
    spdlog::trace("Zoom and Rotate failed.");
    return;
  }

  std::array<Vec3d, 3> relocations;
  std::tie(relocations[0], relocations[2]) = great_prism.value();

  std::vector<std::set<int>> checker;
  auto flag = 1;
  double alpha = 1.;
  while (flag == 1) { // shrink if there is a volume failure.
    Vec3d new_b = relocations[0] * (alpha) + (1 - alpha) * pc.mid[vid];
    Vec3d new_t = relocations[2] * (alpha) + (1 - alpha) * pc.mid[vid];
    flag = prism::local_validity::attempt_relocate(
        pc.base, pc.mid, pc.top, pc.F, *pc.ref.aabb, pc.ref.V, pc.ref.F,
        pc.track_ref, distortion_bound, VF[vid], vid,
        {new_b, pc.mid[vid], new_t}, checker);
    spdlog::trace("newb {}", new_b);
    spdlog::trace("newt {}", new_t);
    alpha *= 0.8;
    if (alpha < 1e-2)
      break;
  }
  if (flag > 0) {
    spdlog::trace("ZoomRotate checker failed.");
    return;
  }
  if (pc.top_grid != nullptr) {
    spdlog::trace("HashGrid remove");
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
                   double distortion_bound, const std::vector<bool> &skip) {
  if (skip[vid])
    return;

  spdlog::trace("smooth attempt: {}", vid);
  auto new_direction = prism::smoother_direction(
      pc.base, pc.mid, pc.top, pc.F, pc.ref.aabb->num_freeze, VF, VFi, vid);

  if (!new_direction) {
    spdlog::trace("No better location.");
    return;
  }

  std::array<Vec3d, 3> relocations{pc.base[vid] + new_direction.value(),
                                   pc.mid[vid] + new_direction.value(),
                                   pc.top[vid] + new_direction.value()};
  auto query =
      [&pc](const Vec3d &s, const Vec3d &t,
            const std::set<int> &total_trackee) -> std::optional<Vec3d> {
    if (pc.ref.aabb->enabled)
      return pc.ref.aabb->segment_query(
          s, t); // this can be discarded if no performance benefit is found.
    std::array<Vec3d, 2> seg_query{s, t};
    for (auto f : total_trackee) {
      auto v0 = pc.ref.F(f, 0), v1 = pc.ref.F(f, 1), v2 = pc.ref.F(f, 2);
      auto mid_intersect = prism::cgal::segment_triangle_intersection(
          seg_query, {pc.ref.V.row(v0), pc.ref.V.row(v1), pc.ref.V.row(v2)});
      if (mid_intersect)
        return mid_intersect;
    }
    return {};
  };

  if (true) { // project onto reference
    std::set<int> total_trackee;
    for (auto f : VF[vid])
      total_trackee.insert(pc.track_ref[f].begin(), pc.track_ref[f].end());
    std::optional<Vec3d> mid_intersect;
    for (int i = 0; i < 20; i++) {
      mid_intersect = query(relocations[0], relocations[2], total_trackee);
      if (mid_intersect)
        break;
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
      pc.base, pc.mid, pc.top, pc.F, *pc.ref.aabb, pc.ref.V, pc.ref.F,
      pc.track_ref, distortion_bound, VF[vid], vid, relocations, checker);

  if (flag > 0) {
    spdlog::trace("Pan checker failed.");
  } else {
    auto &new_tracks = checker;
    for (int i = 0; i < VF[vid].size(); i++) {
      pc.track_ref[VF[vid][i]] = std::move(new_tracks[i]);
    }
    spdlog::trace("SUCCEED, move to next.");
  }
};
} // namespace prism::local
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
  for (int i = 0; i < pc.feature_edges.rows(); i++)
    for (auto j : {0, 1}) {
      if (pc.feature_edges(i, j) < 0)
        continue;
      skip_flag[pc.feature_edges(i, j)] = true;
    }
  for (int i = 0; i < pc.ref.aabb->num_freeze; i++)
    skip_flag[i] = true;
  {
    RowMati mF, mE;
    vec2eigen(pc.F, mF);
    igl::vertex_triangle_adjacency(pc.mid.size(), mF, VF, VFi);
    prism::local::vertex_coloring(mF, groups);
    igl::boundary_facets(mF, mE);
    spdlog::info("boundary, mE {}", mE.rows());
    for (int i = 0; i < mE.rows(); i++)
      for (auto j : {0, 1})
        skip_flag[mE(i, j)] = true;
  }
  spdlog::info("Group Count {}", groups.size());

  std::srand(0);
  if (option.parallel && pc.top_grid!=nullptr)  {
    spdlog::error("Multithread hashmap is not safe. Todo: move to TBB.");
  }
  for (auto &gr : groups)
    igl::parallel_for(gr.size(),
                      [&gr, &pc, &VF, &VFi, &option, &skip_flag](auto ii) {
                        smooth_single(pc, gr[ii], VF, VFi,
                                      option.distortion_bound, skip_flag);
                      },
                      size_t(option.parallel ? 1 : pc.mid.size()));
  for (auto &gr : groups)
    igl::parallel_for(
        gr.size(),
        [&gr, &pc, &VF, &VFi, &option, &skip_flag](auto ii) {
          smooth_prism(pc, gr[ii], VF, VFi, option.distortion_bound,
                       option.target_thickness, false, skip_flag);
          smooth_prism(pc, gr[ii], VF, VFi, option.distortion_bound,
                       option.target_thickness, true, skip_flag);
        },
        size_t(option.parallel ? 1 : pc.mid.size()));

  return;
}

namespace prism::local {

void legacy_smooth_prism(PrismCage &pc, int vid,
                         const std::vector<std::vector<int>> &VF,
                         const std::vector<std::vector<int>> &VFi,
                         double distortion_bound, double target_thick,
                         const std::vector<bool> &skip, bool on_base) {
  if (skip[vid])
    return;

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
      pc.base, pc.mid, pc.top, pc.F, *pc.ref.aabb, pc.ref.V, pc.ref.F,
      pc.track_ref, distortion_bound, VF[vid], vid, relocations, checker);
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
    for (int i = 0; i < pc.ref.aabb->num_freeze; i++)
      skip_flag[i] = true;
  }
  spdlog::info("Group Count {}", groups.size());
  spdlog::set_level(spdlog::level::debug);
  std::srand(0);

  for (auto &gr : groups)
    igl::parallel_for(
        gr.size(),
        [&gr, &pc, &VF, &VFi, &option, &skip_flag](auto ii) {
          legacy_smooth_prism(pc, gr[ii], VF, VFi, option.distortion_bound,
                              option.target_thickness, skip_flag, true);
          legacy_smooth_prism(pc, gr[ii], VF, VFi, option.distortion_bound,
                              option.target_thickness, skip_flag, false);
        },
        size_t(option.parallel ? 1 : pc.mid.size()));

  spdlog::set_level(spdlog::level::info);
  return;
}

} // namespace prism::local
