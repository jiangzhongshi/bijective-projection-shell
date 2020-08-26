#include "remesh_pass.hpp"

#include <igl/boundary_facets.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <queue>

#include "local_mesh_edit.hpp"
#include "prism/PrismCage.hpp"
#include "prism/cage_utils.hpp"
#include "prism/energy/prism_quality.hpp"
#include "prism/feature_utils.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/spatial-hash/AABB_hash.hpp"
#include "retain_triangle_adjacency.hpp"
#include "validity_checks.hpp"

namespace prism::local_validity {
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
    int orig_fnum);
}
namespace collapse {
bool satisfy_link_condition(const std::vector<Vec3i> &F,
                            const std::vector<Vec3i> &FF,
                            const std::vector<Vec3i> &FFi, int f, int e,
                            std::vector<std::pair<int, int>> &neighbor0,
                            std::vector<std::pair<int, int>> &neighbor1) {
  auto f0 = f, e0 = e;
  auto f1 = FF[f0][e0], e1 = FFi[f0][e0];
  auto u0 = F[f0][e0], u1 = F[f1][e1];
  if (f1 == -1)
    return false; // no boundary here
  assert(f1 != -1);
  assert(e1 != -1);
  assert(F[f1][(e1 + 1) % 3] == u0);
  // clockwise
  auto flag0 = prism::get_star_edges(F, FF, FFi, f0, e0, neighbor0);
  auto flag1 = prism::get_star_edges(F, FF, FFi, f1, e1, neighbor1);
  if (!flag0)
    return false;
  std::vector<int> nv0(neighbor0.size()), nv1(neighbor1.size());
  for (int i = 0; i < neighbor0.size(); i++) {
    auto [f, e] = neighbor0[i];
    nv0[i] = F[f][(e + 1) % 3];
  }
  for (int i = 0; i < neighbor1.size(); i++) {
    auto [f, e] = neighbor1[i];
    nv1[i] = F[f][(e + 1) % 3];
  }

  std::sort(nv0.begin(), nv0.end());
  std::sort(nv1.begin(), nv1.end());
  decltype(nv0) inter;
  std::set_intersection(nv0.begin(), nv0.end(), nv1.begin(), nv1.end(),
                        std::back_inserter(inter));
  if (inter.size() == 2)
    return true;
  else
    return false;
}

} // namespace collapse

int prism::local::wildcollapse_pass(PrismCage &pc, RemeshOptions &option) {
  spdlog::enable_backtrace(40);
  auto &F = pc.F;
  auto &V = pc.mid;
  using queue_entry =
      std::tuple<double /*should negative*/, int /*f*/, int /*e*/, int /*u0*/,
                 int /*u1*/, int /*ts*/>;
  std::priority_queue<queue_entry> queue;

  // build connectivity
  auto [FF, FFi] = prism::local::triangle_triangle_adjacency(F);

  // enqueue
  for (auto f = 0; f < F.size(); f++) {
    for (auto e : {0, 1, 2}) {
      auto v0 = F[f][e], v1 = F[f][(e + 1) % 3];
      if (v0 >= v1 || FF[f][e] == -1)
        continue;
      queue.push({-(V[v0] - V[v1]).norm(), f, e, v0, v1, 0});
    }
  }
  std::vector<bool> skip_flag(pc.mid.size(), false);
  for (int i = 0; i < pc.feature_edges.rows(); i++)
    for (auto j : {0, 1}) {
      skip_flag[pc.feature_edges(i, j)] = true;
    }
  for (int i = 0; i < pc.ref.aabb->num_freeze; i++)
    skip_flag[i] = true;

  std::array<int, 5> rejections_steps{0, 0, 0, 0, 0};
  // pop
  int global_tick = 0;
  RowMati timestamp = RowMati::Zero(F.size(), 3);
  while (!queue.empty()) {
    auto [l, f, e, v0, v1, tick] = queue.top();
    l = std::abs(l);
    queue.pop();
    if (f == -1 || FF[f][e] == -1)
      continue;

    auto u0 = F[f][e], u1 = F[f][(e + 1) % 3];
    if (tick != timestamp(f, e))
      continue;
    assert((V[u1] - V[u0]).norm() == l &&
           "Outdated entries will be ignored, this condition can actually "
           "replace the previous");

    if (l * 2.5 > option.sizing_field(V[u0]) * option.target_adjustment[u0] +
                      option.sizing_field(V[u1]) * option.target_adjustment[u1])
      continue; // skip if l > 4/5*(s1+s2)/2

    spdlog::trace("LinkCondition check {} {}", f, e);
    // collapse and misc checks.
    std::vector<std::pair<int, int>> n0, n1;
    if (!collapse::satisfy_link_condition(F, FF, FFi, f, e, n0, n1)) {
      rejections_steps[0]++;
      continue;
    }
    spdlog::debug("n0 {} n1 {}", n0, n1);
    spdlog::trace("LinkCondition pass, attempt {} {}", f, e);
    auto f1 = FF[f][e], e1 = FFi[f][e];
    decltype(n0) n0_, n1_;
    assert(collapse::satisfy_link_condition(F, FF, FFi, f1, e1, n0_, n1_));
    spdlog::debug("_n0 {} n1 {}", n0_, n1_);
    std::tuple<std::vector<int> /*newfid*/, std::vector<int> /*shifts*/,
               std::vector<std::set<int>> /*track*/>
        checker;
    auto repeat_collapse_with_shrinking = [&pc, &option, &skip_flag, &V, &F,
                                           &n0, &n1, f = f, &f1, &u0, &u1,
                                           &checker]() {
      auto flag = 1;
      std::tuple<Vec3d, Vec3d> recover_coordinates{pc.base[u1], pc.top[u1]};
      for (auto rp = 0; rp < 5; rp++) {
        flag = skip_flag[u0]
                   ? 4
                   : prism::local_validity::attempt_collapse(
                         pc.base, V, pc.top, F, *pc.ref.aabb,
                         {pc.base_grid, pc.top_grid}, pc.ref.V, pc.ref.F,
                         pc.track_ref, option.distortion_bound,
                         option.collapse_quality_threshold, n0, n1, f, f1, u0,
                         u1, checker);
        if (flag != 1)
          break;
        pc.base[u1] = (pc.base[u1] + pc.mid[u1]) / 2;
        pc.top[u1] = (pc.top[u1] + pc.mid[u1]) / 2;
      }
      if (flag != 0)
        std::tie(pc.base[u1], pc.top[u1]) = recover_coordinates;
      return flag;
    };
    auto flag = repeat_collapse_with_shrinking();

    spdlog::trace("Attempt Collapse, {} {} pass: {}", f, e, flag == 0);
    if (flag != 0) {
      rejections_steps[flag]++;
      // test the reverse
      std::swap(n0, n1);
      std::swap(f, f1);
      std::swap(u0, u1);
      std::swap(v0, v1);
      std::swap(e, e1);
      decltype(n0) n0_, n1_;
      assert(collapse::satisfy_link_condition(
          F, FF, FFi, f, e, n0_, n1_)); // the previous is already tested
      flag = repeat_collapse_with_shrinking();
      if (flag != 0)
        continue; // still failing
    }
    auto &[new_fid, new_shifts, new_tracks] = checker;
    assert(new_fid.size() == new_shifts.size());
    prism::edge_collapse(F, FF, FFi, f, e);
    spdlog::trace("EdgeCollapse done {} {}", f, e);

    if (pc.top_grid != nullptr) {
      spdlog::trace("HashGrid Update");
      for (auto [f, e] : n0) {
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
    // Not removing replaced ones since (v0,v1) will serve as timestamp.

    global_tick++;
    for (int i = 0; i < new_fid.size(); i++) {
      auto f = new_fid[i];
      for (auto e = 0; e < 3; e++) {
        auto u0 = F[f][e], u1 = F[f][(e + 1) % 3];
        if (u0 >= u1 || FF[f][e] == -1) {
          timestamp(f, e) = -1;
          continue;
        }
        queue.push({-(V[u1] - V[u0]).norm(), f, e, u0, u1, global_tick});
        timestamp(f, e) = global_tick;
        spdlog::trace("pq {} {} {} {} {}", f, e, u0, u1, global_tick);
      }
    }
    spdlog::trace("Edge Collapsed {} {}", f, e);
  }
  spdlog::info("Pass Collapse total {}. lk{}, v{} i{} d{} q{}", global_tick,
               rejections_steps[0], rejections_steps[1], rejections_steps[2],
               rejections_steps[3], rejections_steps[4]);
  Eigen::VectorXi vid_ind, vid_map;
  pc.cleanup_empty_faces(vid_map, vid_ind);
  for (int i = 0; i < vid_ind.size(); i++) {
    option.target_adjustment[i] = option.target_adjustment[vid_ind[i]];
  }
  option.target_adjustment.resize(vid_ind.size());
  std::map<std::pair<int, int>, std::pair<int, std::vector<int>>> new_metas;
  for (auto m : pc.meta_edges) {
    auto [u0, u1] = m.first;
    new_metas[{vid_map[u0], vid_map[u1]}] = m.second;
  }
  pc.meta_edges = std::move(new_metas);
  return global_tick;
}

namespace prism::local {
int featurecollapse_pass(PrismCage &pc, RemeshOptions &option,
                         std::vector<std::set<int>> &feature_region_segments) {
  auto &meta_edges = pc.meta_edges;
  constexpr auto project_vertex_to_segment =
      [](const RowMatd &refV, int list_size, Vec3d &e0, Vec3d &e1,
         std::vector<Vec3d> &proj) -> bool {
    proj.reserve(list_size);
    for (int i = 0, n = list_size; i <= n - 1; i++) {
      proj.emplace_back((i) * (e1 - e0) / (n - 1) + e0);
    }
    return true; // success flag for snap.
  };
  // meta_edges maps a single edge on the middle surface to a chain of edges on
  // reference.
  auto &F = pc.F;
  auto &V = pc.mid;
  using queue_entry =
      std::tuple<double /*should negative*/, int /*f*/, int /*e*/, int /*u0*/,
                 int /*u1*/, int /*ts*/>;
  std::priority_queue<queue_entry> queue;
  int orig_fnum = F.size();
  int inf_node = pc.mid.size();
  for (int f = 0; f < orig_fnum; f++) {
    for (auto e : {0, 1, 2}) {
      auto v0 = F[f][e], v1 = F[f][(e + 1) % 3];
      if (meta_edges.find({v0, v1}) == meta_edges.end())
        continue;
      queue.push({-(V[v0] - V[v1]).norm(), f, e, v0, v1, 0});
    }
  }
  std::set<int> feature_verts;
  for (auto [k, ignore] : meta_edges) {
    feature_verts.insert(k.first);
    feature_verts.insert(k.second);
  }

  auto [FF, FFi] = prism::local::triangle_triangle_adjacency(F);

  std::array<int, 5> rejections_steps{0, 0, 0, 0, 0};
  // pop
  int global_tick = 0;
  while (!queue.empty()) {
    auto [l, f, e, v0, v1, ignore] = queue.top();
    l = std::abs(l);
    queue.pop();
    if (f == -1 || FF[f][e] == -1)
      continue; // skip collapsed

    auto u0 = F[f][e], u1 = F[f][(e + 1) % 3];
    if (u0 == u1 || u0 != v0 ||
        u1 != v1) // vid changed, means the edge is outdated.
      continue;
    assert((V[u1] - V[u0]).norm() == l &&
           "Outdated entries will be ignored, this condition can actually "
           "replace the previous");

    if (l * 2.5 > option.sizing_field(V[u0]) * option.target_adjustment[u0] +
                      option.sizing_field(V[u1]) * option.target_adjustment[u1])
      continue; // skip if l > 4/5*(s1+s2)/2

    spdlog::trace(">>>>>>LinkCondition check {} {}", f, e);
    // collapse and misc checks.
    std::vector<std::pair<int, int>> n0, n1;
    if (!collapse::satisfy_link_condition(F, FF, FFi, f, e, n0, n1)) {
      rejections_steps[0]++;
      continue;
    }

    spdlog::trace("LinkCondition pass, attempt {} {}", f, e);

    auto f1 = FF[f][e], e1 = FFi[f][e]; // fake face and edge

    // snapping location
    auto seg1 = meta_edges.find({u0, u1});
    assert(seg1 != meta_edges.end() && (seg1->second.second[0]) != -1);
    auto u2 = -1;
    for (auto [f, e] :
         n0) { // find the adjacent (in-chain) vertex of u0 except u1.
      auto v = F[f][(e + 1) % 3];
      if (v != u1 && feature_verts.find(v) != feature_verts.end() &&
          meta_edges.find({v, u0}) != meta_edges.end())
        u2 = v;
    }
    if (u2 == -1) {
      spdlog::trace(">>>> Hit Head ");
      rejections_steps[0]++;
      continue;
    }
    spdlog::trace("u1 {} u0 {} u2 {}", u1, u0, u2);
    auto seg0 = meta_edges.find({u2, u0});
    assert(seg0 != meta_edges.end());

    int chain_id = seg1->second.first;
    spdlog::trace("chain {} ", chain_id);
    if (chain_id != seg0->second.first) {
      spdlog::trace("different chain");
      continue;
    }

    std::vector<int> newseg = seg0->second.second;
    assert(newseg.back() == seg1->second.second.front());
    newseg.pop_back();
    newseg.insert(newseg.end(), seg1->second.second.begin(),
                  seg1->second.second.end());
    std::vector<Vec3d> seg_snaps;
    project_vertex_to_segment(pc.ref.V, newseg.size(), V[u2], V[u1], seg_snaps);

    std::vector<Vec3d> seg_restore(newseg.size());
    for (int i = 0; i < newseg.size(); i++) {
      seg_restore[i] = pc.ref.V.row(newseg[i]);
      pc.ref.V.row(newseg[i]) = seg_snaps[i];
    }
    std::tuple<std::vector<int> /*newfid*/, std::vector<int> /*shifts*/,
               std::vector<std::set<int>> /*track*/>
        checker;
    int flag = prism::local_validity::attempt_boundary_collapse(
        pc.base, V, pc.top, F, *pc.ref.aabb, pc.ref.V, pc.ref.F, pc.track_ref,
        option.distortion_bound, option.collapse_quality_threshold, n0, f, f1,
        u0, u1, u2, chain_id, feature_region_segments, checker, orig_fnum);
    spdlog::trace("Attempt Boundary Collapse, {} {} pass: {}", f, e,
                  (flag == 0));
    if (flag != 0) {
      std::for_each(newseg.begin(), newseg.end(),
                    [](auto &a) { assert(a < 2000); });
      for (auto i = 0; i < newseg.size(); i++)
        pc.ref.V.row(newseg[i]) = seg_restore[i];
      rejections_steps[flag]++;
      continue;
    }

    auto &[new_fid, new_shifts, new_tracks] = checker;
    assert(new_fid.size() == new_shifts.size());

    prism::edge_collapse(F, FF, FFi, f, e);
    spdlog::trace("EdgeCollapse done {} {}", f, e);

    if (pc.top_grid != nullptr) {
      spdlog::trace("HashGrid remove");
      for (auto [f, e] : n0) {
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
    meta_edges[{u2, u1}] = {chain_id, newseg};
    spdlog::debug("insert {} {}", u2, u1);
    meta_edges.erase(seg0);
    meta_edges.erase(seg1);
    for (int i = 0; i < new_fid.size(); i++) {
      auto f = new_fid[i];
      for (auto e = 0; e < 3; e++) {
        auto v0 = F[f][e], v1 = F[f][(e + 1) % 3];
        if (v0 == u2 && v1 == u1)
          queue.push({-(V[v0] - V[v1]).norm(), f, e, v0, v1, global_tick});
      }
    }

    global_tick++;
    spdlog::trace("Edge Collapsed {} {}", f, e);
  }
  spdlog::info("Pass Collapse total {}. lk{}, v{} i{} d{} q{}", global_tick,
               rejections_steps[0], rejections_steps[1], rejections_steps[2],
               rejections_steps[3], rejections_steps[4]);
  F.resize(orig_fnum);
  pc.feature_edges.resize(meta_edges.size(), 2);
  int i = 0;
  for (auto [m, e] : meta_edges) {
    pc.feature_edges(i, 0) = m.first;
    pc.feature_edges(i, 1) = m.second;
    i++;
  }
  Eigen::VectorXi vid_map, vid_ind; // new to old
  pc.cleanup_empty_faces(vid_map, vid_ind);
  for (int i = 0; i < vid_ind.size(); i++) {
    option.target_adjustment[i] = option.target_adjustment[vid_ind[i]];
  }
  option.target_adjustment.resize(vid_ind.size());
  std::map<std::pair<int, int>, std::pair<int, std::vector<int>>> new_metas;
  for (auto m : meta_edges) {
    auto [u0, u1] = m.first;
    new_metas[{vid_map[u0], vid_map[u1]}] = m.second;
  }
  meta_edges = std::move(new_metas);
  return global_tick;
}
} // namespace prism::local