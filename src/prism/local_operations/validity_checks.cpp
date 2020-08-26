#include "validity_checks.hpp"

#include <igl/Timer.h>
#include <igl/volume.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include "local_mesh_edit.hpp"
#include <limits>
#include <prism/predicates/inside_octahedron.hpp>
#include <prism/predicates/triangle_triangle_intersection.hpp>
#include <vector>

#include "../geogram/AABB.hpp"
#include "../energy/map_distortion.hpp"
#include "../energy/prism_quality.hpp"
#include "../predicates/inside_prism_tetra.hpp"
#include "../predicates/positive_prism_volume_12.hpp"
#include "../spatial-hash/AABB_hash.hpp"
#include "local_mesh_edit.hpp"

bool dynamic_intersect_check(const std::vector<Vec3d> &base,
                             const std::vector<Vec3i> &F,
                             const std::vector<int> &removed,
                             const std::vector<Vec3i> &tris,
                             const prism::HashGrid &grid) {
  spdlog::debug("In DIC 2x{}", tris.size());
  constexpr auto share_vertex = [](const auto &F, int c, auto &f) {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (F[c][i] == f[j])
          return true;
    return false; // TODO this is not complete, one vertex touch requires
                  // reduced collision check.
  };
  constexpr auto intersection = [](const auto &V, const auto &F, int c,
                                   auto &f) {
    auto &[p0, p1, p2] = F[c];
    auto &[q0, q1, q2] = f;
    return prism::predicates::triangle_triangle_overlap(
        {V[p0], V[p1], V[p2]}, {V[q0], V[q1], V[q2]});
  };
  igl::Timer timer;
  timer.start();
  for (auto f : tris) {
    RowMat3d local;
    for (auto k : {0, 1, 2})
      local.row(k) = base[f[k]];
    std::set<int> candidates;
    grid.query(local.colwise().minCoeff(), local.colwise().minCoeff(),
               candidates);
    std::vector<int> result;
    std::set_difference(candidates.begin(), candidates.end(), removed.begin(),
                        removed.end(), std::inserter(result, result.end()));
    for (auto c : result) {
      if ((!share_vertex(F, c, f)) && intersection(base, F, c, f))
        return false;
    }
  }
  auto elapsed = timer.getElapsedTimeInMicroSec();
  spdlog::debug("DIC true {}", elapsed);
  return true; // safe operation, no intersection
}

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
    if (std::isnan(q))
      return std::numeric_limits<double>::infinity();
    quality = std::max(quality, q);
  }
  return quality;
}

constexpr auto triangle_shifts = [](auto &moved_tris) {
  std::vector<int> new_shifts(moved_tris.size());
  for (int i = 0; i < moved_tris.size(); i++) {
    auto [s, mt, shift] = tetra_split_AorB(moved_tris[i]);
    moved_tris[i] = mt;
    new_shifts[i] = shift;
  }
  return std::move(new_shifts);
};

bool prism_positivity_with_numerical(const std::array<Vec3d, 6> &verts,
                                     const std::array<bool, 3> &constrained) {
  auto matV = Eigen::Map<const RowMatd>(verts[0].data(), 6, 3);
  for (int p = 0; p < 3; p++) {
    if (constrained[p])
      continue;
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
    const std::vector<Vec3d> &mid, // placed new verts
    const std::vector<Vec3d> &top, const std::vector<Vec3i> &tris,
    const std::set<int> &combined_trackee, // indices to ref.F tracked
    const RowMatd &refV, const RowMati &refF, double distortion_bound,
    int num_freeze, bool bundled_intersection = false) {
  spdlog::trace("In DC ct#{}, tris{}", combined_trackee.size(), tris.size());
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
    std::array<bool, 3> oct_type;
    prism::determine_convex_octahedron(base_vert, top_vert, oct_type,
                                       num_freeze > v0);
    spdlog::trace("checking tris{}: {}-{}-{}", i, v0, v1, v2);
    for (auto t : combined_trackee) { // for every tracked original triangle.
      auto ref_tri = std::array<Vec3d, 3>{
          refV.row(refF(t, 0)), refV.row(refF(t, 1)), refV.row(refF(t, 2))};
      bool intersected_prism = false;
      if (v0 >= num_freeze || v0 != refF(t, 0))
        intersected_prism =
            prism::triangle_intersect_octahedron(base_vert, mid_vert, oct_type,
                                                 ref_tri, num_freeze > v0) ||
            prism::triangle_intersect_octahedron(mid_vert, top_vert, oct_type,
                                                 ref_tri, num_freeze > v0);
      else
        intersected_prism = prism::singularless_triangle_intersect_octahedron(
                                base_vert, mid_vert, oct_type, ref_tri) ||
                            prism::singularless_triangle_intersect_octahedron(
                                mid_vert, top_vert, oct_type, ref_tri);
      spdlog::trace("checking vs ref{}: {}", t, intersected_prism);
      if (!intersected_prism)
        continue;

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
      if (bundled_intersection && // this is enabled when global (AABB tree)
                                  // intersection is off.
          (prism::predicates::triangle_triangle_overlap(base_vert,
                                                             ref_tri) ||
           prism::predicates::triangle_triangle_overlap(top_vert,
                                                             ref_tri))) {
        return {};
      }
    }
  }

  auto elapsed = timer.getElapsedTimeInMicroSec();
  spdlog::trace("DC true {}", elapsed);

  //
  return distributed_refs;
}

// below is the checks for local operations

int attempt_collapse(
    const std::vector<Vec3d> &base, const std::vector<Vec3d> &mid,
    const std::vector<Vec3d> &top, const std::vector<Vec3i> &F,
    const prism::geogram::AABB &tree,
    const std::array<std::shared_ptr<prism::HashGrid>, 2> &grid,
    const RowMatd &refV, const RowMati &refF,
    const std::vector<std::set<int>> &map_track, double distortion_bound,
    double improve_quality_threshold,
    // specified infos below
    const std::vector<std::pair<int, int>> &neighbor0,
    const std::vector<std::pair<int, int>> &neighbor1, int f0, int f1,
    int u0, int u1,
    std::tuple<std::vector<int> /*newfid*/, std::vector<int> /*shifts*/,
               std::vector<std::set<int>> /*track*/> &checker) {
  std::vector<Vec3i> moved_tris, old_tris;
  std::vector<int> new_fid, old_fid;
  moved_tris.reserve(neighbor0.size() - 2);
  for (auto [f, e] : neighbor0) {
    assert(F[f][e] == u0);
    auto new_tris = F[f];
    old_fid.push_back(f);
    if (new_tris[0] == u1 || new_tris[1] == u1 || new_tris[2] == u1)
      continue; // collapsed faces
    old_tris.push_back(new_tris);
    new_tris[e] = u1;
    moved_tris.emplace_back(new_tris);
    new_fid.push_back(f);
  }
  // note that new_fid is geometric cyclic.
  assert(moved_tris.size() == neighbor0.size() - 2);

  spdlog::trace("Quality check");
  auto quality_before = max_quality_on_tris(base, mid, top, old_tris);
  auto quality_after = max_quality_on_tris(base, mid, top, moved_tris);
  spdlog::trace("Quality compare {} -> {}", quality_before, quality_after);
  if (std::isnan(quality_after))
    return 4;
  if ((quality_after > improve_quality_threshold) &&
      quality_after > quality_before) // if the quality is too large, not allow
                                      // it to increase.
    return 4;
  //  volume check
  auto vc = volume_check(base, mid, top, moved_tris, tree.num_freeze);
  if (!vc)
    return 1;
  bool ic = false;
  auto [base_grid, top_grid] = grid;
  if (top_grid == nullptr)
    ic = intersect_check(base, top, moved_tris, tree);
  else {
    ic = dynamic_intersect_check(base, F, old_fid, moved_tris, *base_grid) &&
         dynamic_intersect_check(top, F, old_fid, moved_tris, *top_grid);
  }
  if (!ic)
    return 2;
  // get new subdivide types
  auto new_shifts = triangle_shifts(moved_tris);

  std::set<int> combined_tracks;
  for (auto [f, e] : neighbor0) {
    std::merge(map_track[f].begin(), map_track[f].end(),
               combined_tracks.begin(), combined_tracks.end(),
               std::inserter(combined_tracks, combined_tracks.begin()));
  }
  auto sub_refs = distort_check(base, mid, top, moved_tris, combined_tracks,
                                refV, refF, distortion_bound, tree.num_freeze,
                                ! tree.enabled);
  if (!sub_refs)
    return 3;

  checker = std::tuple(std::move(new_fid), std::move(new_shifts),
                       std::move(sub_refs.value()));
  return 0;
}

int attempt_flip(const std::vector<Vec3d> &base, const std::vector<Vec3d> &mid,
                 const std::vector<Vec3d> &top, const std::vector<Vec3i> &F,
                 const prism::geogram::AABB &tree, const RowMatd &refV,
                 const RowMati &refF,
                 const std::vector<std::set<int>> &map_track,
                 double distortion_bound,
                 // specified infos below
                 int f0, int f1, int e0, int e1, int v0, int v1,
                 std::tuple<std::vector<int>,             /*shift*/
                            std::vector<std::set<int>> /*track*/
                            > &checker) {
  std::vector<Vec3i> moved_tris{F[f0], F[f1]};

  auto quality_before = max_quality_on_tris(base, mid, top, moved_tris);
  auto e01 = (e0 + 1) % 3;
  auto e11 = (e1 + 1) % 3;
  moved_tris[0][e01] = v1;
  moved_tris[1][e11] = v0;
  auto quality_after = max_quality_on_tris(base, mid, top, moved_tris);
  spdlog::trace("q{} -> q{}", quality_before, quality_after);
  if (quality_after > quality_before || std::isnan(quality_after))
    return 4;

  //  volume check
  auto vc = volume_check(base, mid, top, moved_tris, tree.num_freeze);
  if (!vc)
    return 1;
  auto ic = intersect_check(base, top, moved_tris, tree);
  if (!ic)
    return 2;

  // get new subdivide types
  auto new_shifts = triangle_shifts(moved_tris);

  std::set<int> combined_tracks;
  std::merge(map_track[f0].begin(), map_track[f0].end(), map_track[f1].begin(),
             map_track[f1].end(),
             std::inserter(combined_tracks, combined_tracks.begin()));

  auto sub_refs = distort_check(base, mid, top, moved_tris, combined_tracks,
                                refV, refF, distortion_bound, tree.num_freeze,! tree.enabled);
  if (!sub_refs)
    return 3;

  checker = std::tuple(std::move(new_shifts), std::move(sub_refs.value()));
  return 0;
}

int attempt_relocate(std::vector<Vec3d> &base, std::vector<Vec3d> &mid,
                     std::vector<Vec3d> &top, const std::vector<Vec3i> &F,
                     const prism::geogram::AABB &tree, const RowMatd &refV,
                     const RowMati &refF,
                     const std::vector<std::set<int>> &map_track,
                     double distortion_bound,
                     // specified infos below
                     const std::vector<int> &nb, int vid,
                     const std::array<Vec3d, 3> /*b-m-t*/ &relocations,
                     std::vector<std::set<int>> &sub_trackee) {
  std::tuple<Vec3d, Vec3d, Vec3d> stored_loc = {base[vid], mid[vid], top[vid]};
  std::vector<Vec3i> nb_tris(nb.size());
  for (int i = 0; i < nb.size(); i++)
    nb_tris[i] = F[nb[i]];
  auto quality_before = max_quality_on_tris(base, mid, top, nb_tris);
  base[vid] = relocations[0];
  mid[vid] = relocations[1];
  top[vid] = relocations[2];

  auto quality_after = max_quality_on_tris(base, mid, top, nb_tris);
  if (quality_after > quality_before || std::isnan(quality_after)) {
    std::tie(base[vid], mid[vid], top[vid]) = stored_loc;
    return 4;
  }

  //  volume check
  if (auto vol_good = volume_check(base, mid, top, nb_tris, tree.num_freeze);
      !vol_good || !intersect_check(base, top, nb_tris, tree)) {
    // assert(vol_good);
    std::tie(base[vid], mid[vid], top[vid]) = stored_loc;
    spdlog::trace("failed inter check {}, with volc={}, restore", vid,
                  vol_good);
    if (!vol_good)
      return 1;
    else
      return 2;
  }

  // get new subdivide types
  auto new_shifts = triangle_shifts(nb_tris);
  for (auto &s : new_shifts)
    assert(s == 0);

  std::set<int> combined_tracks;
  for (auto f : nb) {
    std::merge(map_track[f].begin(), map_track[f].end(),
               combined_tracks.begin(), combined_tracks.end(),
               std::inserter(combined_tracks, combined_tracks.begin()));
  }

  auto sub_refs = distort_check(base, mid, top, nb_tris, combined_tracks, refV,
                                refF, distortion_bound, tree.num_freeze,
                                ! tree.enabled);
  if (!sub_refs) {
    std::tie(base[vid], mid[vid], top[vid]) = stored_loc;
    spdlog::trace("failed map check {}, restore", vid);
    return 3;
  }

  sub_trackee = std::move(sub_refs.value());
  return 0;
}

int attempt_split(
    std::vector<Vec3d> &base, std::vector<Vec3d> &mid, std::vector<Vec3d> &top,
    const std::vector<Vec3i> &F, const prism::geogram::AABB &tree,
    const RowMatd &refV, const RowMati &refF,
    const std::vector<std::set<int>> &map_track, double distortion_bound,
    bool improve_quality,
    // specified infos below
    int f0, int f1, int e0, int e1,
    const std::array<Vec3d, 3> /*b-m-t*/ &newlocations,
    std::tuple<std::vector<int> /*fid*/, std::vector<int>, /*shift*/
               std::vector<std::set<int>>               /*track*/
               > &checker) {
  base.push_back(newlocations[0]);
  mid.push_back(newlocations[1]);
  top.push_back(newlocations[2]);

  int ux = mid.size() - 1;
  std::vector<Vec3i> new_tris{F[f0], F[f1], F[f0], F[f1]};
  // conform to bool edge_split
  new_tris[0][e0] = ux;
  new_tris[1][e1] = ux;
  new_tris[2][(e0 + 1) % 3] = ux;
  new_tris[3][(e1 + 1) % 3] = ux;
  auto quality_before = max_quality_on_tris(base, mid, top, {F[f0], F[1]});
  auto quality_after = max_quality_on_tris(base, mid, top, new_tris);
  spdlog::trace("Quality compare {} -> {}", quality_before, quality_after);
  if (std::isnan(quality_after) ||
      improve_quality && quality_after > quality_before) {
    base.pop_back();
    mid.pop_back();
    top.pop_back();
    return 4;
  }
  // if (quality_after > 1e3 && quality_after > quality_before) return 4;

  std::vector<int> new_fid{f0, f1, static_cast<int>(F.size()),
                           static_cast<int>(F.size() + 1)};

  //  volume check
  auto vc = volume_check(base, mid, top, new_tris, tree.num_freeze);
  if (!vc || !intersect_check(base, top, new_tris, tree)) {
    base.pop_back();
    mid.pop_back();
    top.pop_back();
    if (!vc) {
      spdlog::trace("Split Vol fail");
      return 1;
    } else
      spdlog::trace("Split: Failed inter check with vol={}", vc);
    return 2;
  }

  auto new_shifts = triangle_shifts(new_tris);
  std::vector<std::set<int>> sub_refs;
  std::set<int> combined_tracks;
  std::merge(map_track[f0].begin(), map_track[f0].end(), map_track[f1].begin(),
             map_track[f1].end(),
             std::inserter(combined_tracks, combined_tracks.begin()));

  auto sub_refs_optional =
      distort_check(base, mid, top, new_tris, combined_tracks, refV, refF,
                    distortion_bound, tree.num_freeze,
                    ! tree.enabled);
  if (sub_refs_optional)
    sub_refs = sub_refs_optional.value();
  if (sub_refs.size() == 0) {
    base.pop_back();
    mid.pop_back();
    top.pop_back();
    spdlog::trace("failed map check f{}e{}, restore", f0, e0);
    return 3;
  }

  checker = std::tuple(std::move(new_fid), std::move(new_shifts),
                       std::move(sub_refs));
  return 0;
}

int attempt_boundary_collapse(
    const std::vector<Vec3d> &base, const std::vector<Vec3d> &mid,
    const std::vector<Vec3d> &top, const std::vector<Vec3i> &F,
    const prism::geogram::AABB &tree, RowMatd &refV, const RowMati &refF,
    const std::vector<std::set<int>> &map_track, double distortion_bound,
    double improve_quality_threshold,
    // specified infos below
    const std::vector<std::pair<int, int>> &neighbor0, int f0, int f1,
    int u0, int u1, int u2, int ci,
    std::vector<std::set<int>> &feature_region_segments,
    std::tuple<std::vector<int> /*newfid*/, std::vector<int> /*shifts*/,
               std::vector<std::set<int>> /*track*/> &checker,
    int orig_fnum) {
  std::vector<Vec3i> moved_tris, old_tris;
  std::vector<int> new_fid;
  moved_tris.reserve(neighbor0.size() - 1);
  auto split2 = -1, split1 = -1;
  for (auto [f, e] : neighbor0) {
    if (f >= orig_fnum)
      continue; // skip fake faces.
    assert(F[f][e] == u0);
    auto new_tris = F[f];
    if (F[f][(e + 1) % 3] == u2)
      split2 = new_fid.size();
    if (F[f][(e + 1) % 3] == u1)
      split1 = new_fid.size();

    if (new_tris[0] == u1 || new_tris[1] == u1 || new_tris[2] == u1)
      continue; // there are two triangles to skip
    old_tris.push_back(new_tris);
    new_tris[e] = u1;
    moved_tris.emplace_back(new_tris);
    new_fid.push_back(f);
  }
  // note that new_fid is geometric cyclic.
  /// assuming clockwise.
  std::vector<int> reverse_f, forward_f;
  int nf = new_fid.size();
  for (int j = split2 + 1; j < split1 + (split1 < split2 ? nf : 0); j++)
    forward_f.push_back(j % nf);
  for (int j = split1; j <= split2 + (split1 < split2 ? 0 : nf); j++)
    reverse_f.push_back(j % nf);
  spdlog::trace("Quality check");
  auto quality_before = max_quality_on_tris(base, mid, top, old_tris);
  auto quality_after = max_quality_on_tris(base, mid, top, moved_tris);
  spdlog::trace("Quality compare {} -> {}", quality_before, quality_after);
  if (std::isnan(quality_after))
    return 4;
  if ((quality_after > improve_quality_threshold) &&
      quality_after > quality_before) // if the quality is too large, not allow
                                      // it to increase.
    return 4;
  //  volume check
  auto vc = volume_check(base, mid, top, moved_tris, tree.num_freeze);
  if (!vc) {
    return 1;
  }

  // auto ic = intersect_check(base, top, moved_tris, tree);
  // if (!ic) return 2;
  // get new subdivide types
  auto new_shifts = triangle_shifts(moved_tris);

  std::set<int> combined_trackee, forward_combined, reverse_combined;
  decltype(moved_tris) moved_f, moved_r;
  for (auto i : forward_f) {
    auto f = new_fid[i];
    moved_f.push_back(moved_tris[i]);
  }
  for (auto i : reverse_f) {
    auto f = new_fid[i];
    moved_r.push_back(moved_tris[i]);
  }
  for (auto [f, e] : neighbor0)
    // for (auto f : new_fid)
    std::merge(map_track[f].begin(), map_track[f].end(),
               combined_trackee.begin(), combined_trackee.end(),
               std::inserter(combined_trackee, combined_trackee.begin()));
  std::set_difference(
      combined_trackee.begin(), combined_trackee.end(),
      feature_region_segments[2 * ci + 1].begin(),
      feature_region_segments[2 * ci + 1].end(),
      std::inserter(forward_combined, forward_combined.begin()));
  std::set_difference(
      combined_trackee.begin(), combined_trackee.end(),
      feature_region_segments[2 * ci].begin(),
      feature_region_segments[2 * ci].end(),
      std::inserter(reverse_combined, reverse_combined.begin()));
  assert(forward_combined.size() > 0);
  assert(reverse_combined.size() > 0);
  auto sub_refs_f =
      distort_check(base, mid, top, moved_f, forward_combined, refV, refF,
                    distortion_bound, tree.num_freeze, /*bundle intersect*/true);
  auto sub_refs_r =
      distort_check(base, mid, top, moved_r, reverse_combined, refV, refF,
                    distortion_bound, tree.num_freeze, /*bundle intersect*/true);
  if (!sub_refs_f || !sub_refs_r)
    return 3;

  std::vector<std::set<int>> sub_refs(new_fid.size());
  for (auto i = 0; i < forward_f.size(); i++) {
    sub_refs[forward_f[i]] = std::move(sub_refs_f.value()[i]);
  }
  for (auto i = 0; i < reverse_f.size(); i++) {
    sub_refs[reverse_f[i]] = std::move(sub_refs_r.value()[i]);
  }
#ifndef NDEBUG
  for (auto &s : sub_refs)
    assert(s.size() > 0);
#endif
  checker = std::tuple(std::move(new_fid), std::move(new_shifts),
                       std::move(sub_refs));
  return 0;
}

} // namespace prism::local_validity