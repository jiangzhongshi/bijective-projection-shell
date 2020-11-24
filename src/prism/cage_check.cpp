#include "cage_check.hpp"

#include <igl/triangle_triangle_adjacency.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/ranges.h>

#include "PrismCage.hpp"
#include "geogram/AABB.hpp"
#include "igl/write_triangle_mesh.h"
#include "prism/energy/map_distortion.hpp"
#include "prism/local_operations/remesh_pass.hpp"
#include "prism/local_operations/retain_triangle_adjacency.hpp"
#include "prism/predicates/inside_octahedron.hpp"

namespace prism::local_validity {
std::optional<std::vector<std::set<int>>> distort_check(
    const std::vector<Vec3d>&,
    const std::vector<Vec3d>&,  // placed new verts
    const std::vector<Vec3d>&, const std::vector<Vec3i>& shells,
    const std::set<int>& combined_trackee,  // indices to ref.F tracked
    const RowMatd& refV, const RowMati& refF, double distortion_bound,
    int num_freeze, bool bundled_intersection);
}
bool prism::cage_check::verify_edge_based_track(
    const PrismCage& pc, const prism::local::RemeshOptions& option,
    std::vector<std::set<int>>& map_track) {
  auto num_freeze = pc.ref.aabb->num_freeze;
  assert(map_track.size() == pc.F.size() && "use with shell pipeline");
  auto verify = [& F_sh = std::as_const(pc.F), &top = std::as_const(pc.top),
                 &base = std::as_const(pc.base), &mid = std::as_const(pc.mid),
                 &refV = pc.ref.V, &refF = pc.ref.F, &option,
                 num_freeze](int sh, int i, double distortion) -> bool {
    auto tracks = prism::local_validity::distort_check(
        base, mid, top, {F_sh[sh]}, std::set<int>{i}, refV, refF,
        distortion, num_freeze, option.dynamic_hashgrid);
        return (tracks && tracks.value()[0].size() > 0);
  };

  // do feature
  std::set<std::pair<int, int>> orig_chains;
  for (auto [m, id_chain] : pc.meta_edges) {
    auto [id, chain] = id_chain;
    // auto v0 = m.first, v1 = m.second;
    // if (v0 > v1) std::swap(v0, v1);
    // orig_chains.insert({v0, v1});
    for (int i = 0; i < chain.size() - 1; i++) {
      auto v0 = chain[i], v1 = chain[(i + 1)];
      if (v0 > v1) std::swap(v0, v1);
      orig_chains.insert({v0, v1});
    }
  }

  // verify that existing shell-triangle pair are valid
  // verify that an extra triangle will not be tracked
  RowMati TT, TTi;
  igl::triangle_triangle_adjacency(pc.ref.F, TT, TTi);
  for (int si = 0; si < pc.F.size(); si++) {
    if (pc.F[si][0] == pc.F[si][1]) continue; // uncleaned collapsed faces
    std::set<int> extended, rims;
    spdlog::debug("track {}", map_track[si]);
    for (auto tri : map_track[si]) {
      spdlog::trace("s t {} {}", si, tri);
      if (verify(si, tri, 0.) != true) {
        spdlog::dump_backtrace();
        spdlog::warn("shell {} SHOULD intersect triangle {}", si, tri);
        return false;
      }
      for (int j = 0; j < 3; j++) {
        auto v0 = pc.ref.F(tri,j), v1 = pc.ref.F(tri,(j + 1) % 3);
        if (v0 > v1) std::swap(v0, v1);
        if (orig_chains.find({v0, v1}) != orig_chains.end()) continue;
        if (TT(tri, j) != -1) extended.insert(TT(tri, j));
      }
      set_minus(extended, map_track[si], rims);
      for (auto tri : rims) {
        spdlog::trace("s t {} not {}", si, tri);
        if (verify(si, tri, option.distortion_bound) != false) {
          spdlog::dump_backtrace();
          spdlog::warn("shell {} should NOT intersect triangle {}", si,
                       tri);
          return false;
        }
      }
    }
  }
  return true;
}
bool prism::cage_check::verify_bijection(
    const PrismCage& pc, std::vector<Vec3d>& V, std::vector<Vec3i>& F,
    std::vector<std::set<int>>& track_to_prism) {
  assert(track_to_prism.size() == F.size() && "for use during section mesh");
  // TODO: refactor and use the section distort check directly.
  auto verify = [& F_sh = std::as_const(pc.F), &top = std::as_const(pc.top),
                 &base = std::as_const(pc.base), &mid = std::as_const(pc.mid),
                 &V, &tris = F](int sh, int i, bool require_inter) -> bool {
    auto cur_tri =
        std::array<Vec3d, 3>{V[tris[i][0]], V[tris[i][1]], V[tris[i][2]]};

    int num_freeze = 0;
    auto [v0, v1, v2] = F_sh[sh];
    std::array<Vec3d, 3> base_vert{base[v0], base[v1], base[v2]};
    std::array<Vec3d, 3> mid_vert{mid[v0], mid[v1], mid[v2]};
    std::array<Vec3d, 3> top_vert{top[v0], top[v1], top[v2]};
    std::array<bool, 3> oct_type;
    prism::determine_convex_octahedron(base_vert, top_vert, oct_type,
                                       num_freeze > v0);

    bool intersected_prism =
        prism::triangle_intersect_octahedron(base_vert, mid_vert, oct_type,
                                             cur_tri, num_freeze > v0) ||
        prism::triangle_intersect_octahedron(mid_vert, top_vert, oct_type,
                                             cur_tri, num_freeze > v0);
    if (!intersected_prism) {  // if no intersection
      if (require_inter) {
        spdlog::trace("require intersection s{} t{}", sh, i);
        return false;
      }
      return true;
    }

    if (require_inter == false) {
      spdlog::trace("unintended intersection occur s{} t{}", sh, i);
      return false;
    }
    for (int tc = (v0 < num_freeze) ? 1 : 0; tc < 3; tc++) {
      auto pillar = top_vert[tc] - base_vert[tc];
      auto distortion = prism::energy::map_max_cos_angle(pillar, cur_tri);
      if (distortion < 0.) {
        spdlog::trace("distortion fail {} occur s{} t{}", distortion, sh, i);
        return false;
      }
    }
    return true;
  };
  auto [TT, TTi] = prism::local::triangle_triangle_adjacency(pc.F);
  for (int i = 0; i < track_to_prism.size();
       i++) {  // for each triangle in section mesh
    auto& si = track_to_prism[i];
    for (auto sf : si) {
      // verify shell sf intersects triangle i;

      if (!verify(sf, i, true)) {
        spdlog::trace("verify shell {} SHOULD intersect triangle {}", sf, i);
        return false;
      }
      for (int j = 0; j < 3; j++) {
        auto f1 = TT[sf][j];
        if (si.find(f1) != si.end()) continue;  // skip if neighbor also in.

        if (!verify(f1, i, false)) {
          spdlog::trace("verify shell {} NOT intersect triangle {}", f1, i);
          return false;
        }
      }
    }
  }
  return true;
}

bool prism::cage_check::cage_is_away_from_ref(PrismCage& pc) {
  auto aabb = new prism::geogram::AABB(pc.ref.V, pc.ref.F);
  for (auto [v0, v1, v2] : pc.F) {
    // singular edge does not exist. Bevel always split it aggressively.
    assert(!(v1 < aabb->num_freeze && v2 < aabb->num_freeze));
    if (aabb->intersects_triangle({pc.base[v0], pc.base[v1], pc.base[v2]},
                                  v0 < aabb->num_freeze)) {
      spdlog::dump_backtrace();
      spdlog::error("Intersect Base [{}, {}, {}]", v0, v1, v2);
      return false;
    }
    if (aabb->intersects_triangle({pc.top[v0], pc.top[v1], pc.top[v2]},
                                  v0 < aabb->num_freeze)) {
      spdlog::error("Intersect Top [{}, {}, {}]", v0, v1, v2);
      return false;
    }
  }
  return true;
}