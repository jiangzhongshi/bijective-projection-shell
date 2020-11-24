#include "remesh_pass.hpp"

#include <igl/boundary_facets.h>
#include <igl/doublearea.h>
#include <igl/volume.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <prism/energy/smoother_pillar.hpp>
#include <queue>

#include "local_mesh_edit.hpp"
#include "prism/PrismCage.hpp"
#include "prism/cage_utils.hpp"
#include "prism/energy/prism_quality.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/spatial-hash/AABB_hash.hpp"
#include "retain_triangle_adjacency.hpp"
#include "validity_checks.hpp"

namespace prism::local_validity {

int attempt_split(
    const PrismCage &pc, const std::vector<std::set<int>> &map_track,
    const prism::local::RemeshOptions &option,
    // double distortion_bound, bool improve_quality,
    // specified infos below
    int f0, int f1, int e0, int e1,
    std::tuple<std::vector<int> /*fid*/, std::vector<int>, /*shift*/
               std::vector<std::set<int>>                  /*track*/
               > &checker) {
  auto &base = pc.base, &top = pc.top, &mid = pc.mid;
  auto &F = pc.F;
  auto &refV = pc.ref.V;
  auto &refF = pc.ref.F;
  auto &tree = *pc.ref.aabb;

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
      option.split_improve_quality && quality_after > quality_before) {
    return 4;
  }

  std::vector<int> new_fid{f0, f1, static_cast<int>(F.size()),
                           static_cast<int>(F.size() + 1)};

  //  volume check
  auto vc = volume_check(base, mid, top, new_tris, tree.num_freeze);
  if (!vc || !intersect_check(base, top, new_tris, tree)) {
    if (!vc) {
      spdlog::trace("Split Vol fail");
      return 1;
    } else
      spdlog::trace("Split: Failed inter check with vol={}", vc);
    return 2;
  }

  auto new_shifts = prism::local_validity::triangle_shifts(new_tris);
  std::vector<std::set<int>> sub_refs;
  std::set<int> combined_tracks;
  std::merge(map_track[f0].begin(), map_track[f0].end(), map_track[f1].begin(),
             map_track[f1].end(),
             std::inserter(combined_tracks, combined_tracks.begin()));

  auto sub_refs_optional = distort_check(
      base, mid, top, new_tris, combined_tracks, refV, refF,
      option.distortion_bound, tree.num_freeze, option.dynamic_hashgrid);
  if (sub_refs_optional) sub_refs = sub_refs_optional.value();
  if (sub_refs.size() == 0) {
    spdlog::trace("failed map check f{}e{}, restore", f0, e0);
    return 3;
  }

  checker = std::tuple(std::move(new_fid), std::move(new_shifts),
                       std::move(sub_refs));

  return 0;
}

int attempt_flip(const PrismCage &pc,
                 const std::vector<std::set<int>> &map_track,
                 const prism::local::RemeshOptions &option,
                 // specified infos below
                 int f0, int f1, int e0, int e1, int v0, int v1,
                 std::tuple<std::vector<int>,          /*shift*/
                            std::vector<std::set<int>> /*track*/
                            > &checker) {
  auto &base = pc.base, &top = pc.top, &mid = pc.mid;
  auto &F = pc.F;
  auto &refV = pc.ref.V;
  auto &refF = pc.ref.F;
  auto &tree = *pc.ref.aabb;

  std::vector<Vec3i> moved_tris{F[f0], F[f1]};
  std::vector<int> old_fid{f0, f1};

  auto inter_check = [&pc, &old_fid, &moved_tris, &option]() {
    bool ic = false;
    if (!option.dynamic_hashgrid)
      ic = intersect_check(pc.base, pc.top, moved_tris, *pc.ref.aabb);
    else {
      ic = dynamic_intersect_check(pc.base, pc.F, old_fid, moved_tris,
                                   *pc.base_grid) &&
           dynamic_intersect_check(pc.top, pc.F, old_fid, moved_tris,
                                   *pc.top_grid);
    }
    return ic;
  };

  auto quality_before = max_quality_on_tris(base, mid, top, moved_tris);
  auto e01 = (e0 + 1) % 3;
  auto e11 = (e1 + 1) % 3;
  moved_tris[0][e01] = v1;
  moved_tris[1][e11] = v0;
  auto quality_after = max_quality_on_tris(base, mid, top, moved_tris);
  spdlog::trace("q{} -> q{}", quality_before, quality_after);
  if (quality_after > quality_before || std::isnan(quality_after)) return 4;

  //  volume check
  auto vc = volume_check(base, mid, top, moved_tris, tree.num_freeze);
  if (!vc) return 1;
  if (!inter_check()) return 2;

  // get new subdivide types
  auto new_shifts = prism::local_validity::triangle_shifts(moved_tris);

  std::set<int> combined_tracks;
  std::merge(map_track[f0].begin(), map_track[f0].end(), map_track[f1].begin(),
             map_track[f1].end(),
             std::inserter(combined_tracks, combined_tracks.begin()));

  auto sub_refs = distort_check(base, mid, top, moved_tris, combined_tracks,
                                refV, refF, option.distortion_bound,
                                tree.num_freeze, option.dynamic_hashgrid);
  if (!sub_refs) return 3;

  checker = std::tuple(std::move(new_shifts), std::move(sub_refs.value()));
  return 0;
}
}  // namespace prism::local_validity
namespace prism::local {

void wildflip_pass(PrismCage &pc, const RemeshOptions &option) {
  auto &F = pc.F;
  auto &V = pc.mid;
  using queue_entry =
      std::tuple<double, int /*f*/, int /*e*/, int /*u0*/, int /*u1*/, int>;
  std::priority_queue<queue_entry> queue;

  // build connectivity
  auto [FF, FFi] = prism::local::triangle_triangle_adjacency(F);

  std::set<std::pair<int, int>> skip_edges;
  for (auto [m, ignore] : pc.meta_edges) {
    auto [v0, v1] = m;
    skip_edges.insert({std::min(v0, v1), std::max(v0, v1)});
  }

  // enqueue
  for (auto [f, flags] = std::pair(0, RowMati(RowMati::Zero(F.size(), 3)));
       f < F.size(); f++) {
    for (auto e : {0, 1, 2}) {
      auto v0 = F[f][e], v1 = F[f][(e + 1) % 3];
      if (v0 > v1 || flags(f, e) == 1) continue;
      if (FF[f][e] == -1) continue;
      if (skip_edges.find({v0, v1}) != skip_edges.end()) continue;
      queue.push({(V[v0] - V[v1]).norm(), f, e, v0, v1, 0});
      flags(f, e) = 1;
      flags(FF[f][e], FFi[f][e]) = 1;
    }
  }

  std::vector<int> valence(V.size(), 0);
  for (auto f : F)
    for (int j = 0; j < 3; j++) valence[f[j]]++;

  int global_tick = 0;
  RowMati timestamp = RowMati::Zero(F.size(), 3);
  std::vector<int> rejection_steps(8, 0);
  // pop
  while (!queue.empty()) {
    auto [l, f, e, u0, u1, tick] = queue.top();
    queue.pop();
    if (f == -1 || FF[f][e] == -1) continue;  // skip boundary

    if (auto u0_ = F[f][e], u1_ = F[f][(e + 1) % 3];
        u0_ == u1_ || u0_ != u0 ||
        u1_ != u1)  // vid changed, means the edge is outdated.
      continue;
    else
      assert((V[u1_] - V[u0_]).norm() == l &&
             "Outdated entries will be ignored, this condition can actually "
             "replace the previous");

    auto f1 = FF[f][e], e1 = FFi[f][e];
    assert(f1 != -1);
    auto f0 = f, e0 = e;
    auto v0 = F[f0][(e0 + 2) % 3];
    auto v1 = F[f1][(e1 + 2) % 3];

    // check valency energy.
    constexpr auto valence_energy = [](int i0, int i1, int i2, int i3) {
      double e = 0;
      for (auto i : {i0, i1, i2, i3}) e += (i - 6) * (i - 6);
      return e;
    };
    if (valence_energy(valence[u0], valence[v0], valence[u1], valence[v1]) <
        valence_energy(valence[u0] - 1, valence[v0] + 1, valence[u1] - 1,
                       valence[v1] + 1))
      continue;
    std::tuple<std::vector<int>, std::vector<std::set<int>>> checker;
    int flag = prism::local_validity::attempt_flip(pc, pc.track_ref, option, f0,
                                                   f1, e0, e1, v0, v1, checker);

    std::vector<RowMatd> local_cp;
    std::vector<int> n0_n1{f0, f1};
    std::vector<Vec3i> moved_tris{F[f0], F[f1]};
    if (flag == 0) {
      auto e01 = (e0 + 1) % 3;
      auto e11 = (e1 + 1) % 3;
      moved_tris[0][e01] = v1;
      moved_tris[1][e11] = v0;
      prism::local_validity::triangle_shifts(moved_tris);
      if (option.curve_checker.first.has_value() &&
          !(std::any_cast<std::function<bool(
                const PrismCage &, const std::vector<int> &,
                const std::vector<Vec3i> &, decltype(local_cp) &)>>(
              option.curve_checker.first))(pc, n0_n1, moved_tris, local_cp)) {
        flag = 5;
      }
    }
    if (flag > 0) {
      rejection_steps[flag]++;
      continue;
    }
    auto &[new_shifts, new_tracks] = checker;
    std::vector<int> new_fid{f, f1};

    if (!prism::edge_flip(F, FF, FFi, f, e)) {
      rejection_steps[0]++;
      continue;
    }
    if (option.curve_checker.second.has_value())
      (std::any_cast<
          std::function<void(const std::vector<int> &, const std::vector<int> &,
                             const std::vector<RowMatd> &)>>(
          option.curve_checker.second))(n0_n1, {f0, f1}, local_cp);
    if (pc.top_grid != nullptr) {
      spdlog::trace("HashGrid Update");
      for (auto f : new_fid) {
        pc.top_grid->remove_element(f);
        pc.base_grid->remove_element(f);
      }
      pc.top_grid->insert_triangles(pc.top, F, new_fid);
      pc.base_grid->insert_triangles(pc.base, F, new_fid);
    }

    assert(new_fid.size() == new_tracks.size());
    for (int i = 0; i < new_tracks.size(); i++) {
      pc.track_ref[new_fid[i]] = new_tracks[i];
    }

    // shifts
    shift_left(new_fid, new_shifts, F, FF, FFi);

    // Push the modified edges back in the queue
    // Not removing replaced ones since additional checks are in place: v0,v1
    // and norm==l.
    global_tick++;
    auto push_to_queue = [&queue, &F, &V, &timestamp, &global_tick](auto f,
                                                                    auto v) {
      auto e = -1;
      auto face = F[f];
      for (int i = 0; i < 3; i++)
        if (face[i] == v) e = i;
      if (e == -1) spdlog::error("push queue wrong");
      auto u0 = F[f][e], u1 = F[f][(e + 1) % 3];
      if (u0 > u1) {
        return;
      }
      queue.push({(V[u1] - V[u0]).norm(), f, e, u0, u1, global_tick});
    };
    if (v0 > v1)
      push_to_queue(f1, v1);
    else
      push_to_queue(f1, v0);
    if (v0 < u0) push_to_queue(f0, v0);
    if (u0 < v1) push_to_queue(f0, u0);
    if (v1 < u1) push_to_queue(f1, v1);
    if (u1 < v0) push_to_queue(f1, u1);

    valence[v0]++;
    valence[v1]++;
    valence[u0]--;
    valence[u1]--;
  }
  spdlog::info("Flip {} Done, Rej t{} v{} i{} d{} q{}", global_tick,
               rejection_steps[0], rejection_steps[1], rejection_steps[2],
               rejection_steps[3], rejection_steps[4]);
}

void wildsplit_pass(PrismCage &pc, RemeshOptions &option) {
  auto &F = pc.F;
  auto &V = pc.mid;
  auto input_vnum = V.size();
  using queue_entry =
      std::tuple<double, int /*f*/, int /*e*/, int /*u0*/, int /*u1*/>;
  std::priority_queue<queue_entry> queue;

  // build connectivity
  auto [FF, FFi] = prism::local::triangle_triangle_adjacency(F);

  // enqueue
  for (auto [f, flags] = std::pair(0, RowMati(RowMati::Zero(F.size(), 3)));
       f < F.size(); f++) {
    for (auto e : {0, 1, 2}) {
      auto v0 = F[f][e], v1 = F[f][(e + 1) % 3];
      if (v0 > v1 || flags(f, e) == 1) continue;
      if (FF[f][e] == -1) continue;  // skip boundary
      queue.push({(V[v0] - V[v1]).norm(), f, e, v0, v1});
      flags(f, e) = 1;
      flags(FF[f][e], FFi[f][e]) = 1;
    }
  }

  std::array<int, 6> rejections_steps{0, 0, 0, 0, 0, 0};
  // pop
  while (!queue.empty()) {
    auto [l, f0, e0, u0, u1] = queue.top();
    queue.pop();

    if (f0 == -1 || FF[f0][e0] == -1) continue;  // skip boundary
    if (auto u0_ = F[f0][e0], u1_ = F[f0][(e0 + 1) % 3];
        u0_ == u1_ || u0_ != u0 ||
        u1_ != u1)  // vid changed, means the edge is outdated.
      continue;
    // spdlog::debug("l {} with {} {}", l, option.sizing_field(V[u0]) *
    // option.target_adjustment[u0], option.sizing_field(V[u1]) *
    // option.target_adjustment[u1]);
    if (std::abs(l) * 1.5 <
        (option.sizing_field(V[u0]) * option.target_adjustment[u0] +
         option.sizing_field(V[u1]) * option.target_adjustment[u1])) {
      //  spdlog::debug("skip");
      continue;  // skip if l < 4/3*(s1+s2)/2
    }

    auto f1 = FF[f0][e0], e1 = FFi[f0][e0];
    if (f1 == -1) continue;  // boundary check
    auto v0 = F[f0][(e0 + 2) % 3];
    auto v1 = F[f1][(e1 + 2) % 3];
    std::array<Vec3d, 3> newlocation{(pc.base[u0] + pc.base[u1]) / 2,
                                     Vec3d(0, 0, 0),
                                     (pc.top[u0] + pc.top[u1]) / 2};
    auto new_mid = pc.ref.aabb->segment_query(newlocation[0], newlocation[2]);
    if (!new_mid) {
      spdlog::dump_backtrace();
      spdlog::error("split mid failed");
      spdlog::error("Base {}, {}", pc.base[u0], pc.base[u1]);
      spdlog::error("Top {}, {}", pc.top[u0], pc.top[u1]);
      spdlog::error("New {}, {}", newlocation[0], newlocation[2]);
      continue;
    }
    bool edge_on_ref = false;
    newlocation[1] = new_mid.value();

    spdlog::trace("Attempting: {}-{} {}-{} {}->{} {}-{}", f0, e0, f1, e1, u0,
                  u1, v0, v1);
    std::tuple<std::vector<int> /*fid*/, std::vector<int>, /*shift*/
               std::vector<std::set<int>>                  /*track*/
               >
        checker;

    auto alpha = 1.;
    auto flag = 1;
    while (flag == 1) {  // vc problem
      auto new_b = newlocation[0] * (alpha) + (1 - alpha) * newlocation[1];
      auto new_t = newlocation[2] * (alpha) + (1 - alpha) * newlocation[1];
      pc.base.push_back(newlocation[0]);
      pc.mid.push_back(newlocation[1]);
      pc.top.push_back(newlocation[2]);
      flag = prism::local_validity::attempt_split(pc, pc.track_ref, option, f0,
                                                  f1, e0, e1, checker);
      alpha *= 0.8;
      if (flag != 0) {
        pc.base.pop_back();
        pc.mid.pop_back();
        pc.top.pop_back();
      }
      if (alpha < 1e-2) break;
    }
    std::vector<RowMatd> local_cp;
    if (flag == 0) {
      std::vector<Vec3i> new_tris{F[f0], F[f1], F[f0], F[f1]};
      int ux = pc.mid.size() - 1;
      // conform to bool edge_split
      new_tris[0][e0] = ux;
      new_tris[1][e1] = ux;
      new_tris[2][(e0 + 1) % 3] = ux;
      new_tris[3][(e1 + 1) % 3] = ux;
      prism::local_validity::triangle_shifts(new_tris);

      if (option.curve_checker.first.has_value() &&
          !(std::any_cast<std::function<bool(
                const PrismCage &, const std::vector<int> &,
                const std::vector<Vec3i> &, std::vector<RowMatd> &)>>(
              option.curve_checker.first))(pc, {f0, f1}, new_tris, local_cp)) {
        return;
      }
    }
    if (flag != 0) {
      spdlog::trace("Split Attempt Failed {}-{} {}-{}", f0, e0, f1, e1);
      rejections_steps[flag]++;
      continue;
    }
    auto &[new_fid, new_shifts, new_tracks] = checker;
    assert(new_fid.size() == new_tracks.size());
    prism::edge_split(V.size() - 1, F, FF, FFi, f0, e0);

    if (option.curve_checker.second.has_value())
      (std::any_cast<
          std::function<void(const std::vector<int> &, const std::vector<int> &,
                             const std::vector<RowMatd> &)>>(
          option.curve_checker.second))({f0, f1}, new_fid, local_cp);
    if (pc.top_grid != nullptr) {
      spdlog::trace("HashGrid Update");
      for (auto f : {f0, f1}) {
        pc.top_grid->remove_element(f);
        pc.base_grid->remove_element(f);
      }
      pc.top_grid->insert_triangles(pc.top, F, new_fid);
      pc.base_grid->insert_triangles(pc.base, F, new_fid);
    }

    pc.track_ref.resize(F.size());
    for (int i = 0; i < new_tracks.size(); i++) {
      pc.track_ref[new_fid[i]] = new_tracks[i];
    }

    // shifts
    shift_left(new_fid, new_shifts, F, FF, FFi);
    option.target_adjustment.push_back(
        (option.target_adjustment[u0] + option.target_adjustment[u1]) / 2);

    auto push_to_queue = [&queue, &F, &V, input_vnum](auto f, auto v) {
      auto e = -1;
      auto face = F[f];
      for (int i = 0; i < 3; i++)
        if (face[i] == v) e = i;
      if (e == -1) spdlog::error("push queue wrong");
      auto u0 = F[f][e], u1 = F[f][(e + 1) % 3];
      if (u0 > u1 || u1 >= input_vnum) {
        return;
      }

      queue.push({(V[u1] - V[u0]).norm(), f, e, u0, u1});
      spdlog::trace("pushed {} {} {} {}", f, e, u0, u1);
    };

    auto fx0 = F.size() - 2;
    auto fx1 = F.size() - 1;
    spdlog::trace("Consider {} {} {} {}", u0, u1, v0, v1, fx0, fx1);
    if (v0 < u0) push_to_queue(fx0, v0);
    if (v1 < u1) push_to_queue(fx1, v1);
    if (u1 < v0 && new_shifts[0] != 0) push_to_queue(f0, u1);
    if (u0 < v1 && new_shifts[1] != 0) push_to_queue(f1, u0);
    // push_to_queue(fx0, e0);
    // push_to_queue(f1, (e1 + 2) % 3);
    // push_to_queue(fx1, (e1) % 3);
    // push_to_queue(f0, (e0 + 2) % 3);
  }
  spdlog::info("Split Done, Rejections v{} i{} d{} q{}", rejections_steps[1],
               rejections_steps[2], rejections_steps[3], rejections_steps[4]);

  std::set<int> low_quality_vertices;
  for (auto [v0, v1, v2] : F) {
    if (prism::energy::triangle_quality({V[v0], V[v1], V[v2]}) > 20) {
      low_quality_vertices.insert(v0);
      low_quality_vertices.insert(v1);
      low_quality_vertices.insert(v2);
    }
  }
  for (auto v : low_quality_vertices) option.target_adjustment[v] /= (2 * 1.5);
  for (auto &u : option.target_adjustment)
    u = std::max(std::min(1.5 * u, 1.), 1e-5);

  spdlog::info("Post Split Adjustments: low_quality {}/{}",
               low_quality_vertices.size(), V.size());
}
}  // namespace prism::local
