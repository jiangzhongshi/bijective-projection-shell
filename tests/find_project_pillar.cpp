
#include <doctest.h>

#include <geogram/mesh/mesh_AABB.h>
#include <igl/upsample.h>
#include <spdlog/spdlog.h>
#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_reorder.h>

#include <spdlog/fmt/ostr.h>
#include <igl/barycentric_coordinates.h>
#include <igl/embree/EmbreeIntersector.h>
// for the texture transfer with subdivision.
#include "test_common.hpp"
#include <prism/geogram/geogram_utils.hpp>

Vec3d fiber_hit_ref(const std::vector<Vec3d>& base,
                    const std::vector<Vec3d>& top,
                    const igl::embree::EmbreeIntersector& embree, int v0,
                    int v1, int v2, const Vec3d& facet_bary) {
  auto& rising_facets = v1 > v2 ? PRISM_RISING_FACETS_A : PRISM_RISING_FACETS_B;

  std::array<Vec3d, 6> prism_verts{base[v0], base[v1], base[v2],
                                   top[v0],  top[v1],  top[v2]};

  std::array<Vec3d, 4> fiber_points;
  for (int i = 0; i < 4; i++) {
    fiber_points[i] = Vec3d(0, 0, 0);
    auto& rf = rising_facets[i];
    spdlog::trace("rf {}: {} {} {}", i, rf[0], rf[1], rf[2]);
    for (int j = 0; j < 3; j++)
      fiber_points[i] += prism_verts[rf[j]] * facet_bary[j];
    spdlog::trace("fiber {}: {}", i, fiber_points[i]);
  }
  std::vector<igl::Hit> hits;
  for (int i = 0; i < 3; i++) {
    igl::Hit hit_A;
    auto dir = (fiber_points[i + 1] - fiber_points[i]).cast<float>();
    if (dir.squaredNorm()<1e-10) continue;
    bool hitted = embree.intersectBeam(
        fiber_points[i].cast<float>(), dir, hit_A);
    if (!hitted) continue;
    if (hit_A.t <= 1) {
      spdlog::trace("Hitted: ") ;
      return Vec3d(hit_A.id, hit_A.u, hit_A.v);
    }
    hits.emplace_back(hit_A);
  }
  if (hits.size() == 0) {
    assert("No numerical hit. Blame embree");
    return {};
  }
  else {
    return {};
  }
};

// TEST_CASE("find pillar") {
//   spdlog::set_level(spdlog::level::warn);
//   PrismCage pc("../tests/It2_cfs.h5");
//   RowMatd subV, mV, refV;
//   RowMati subF, mF, refF;
//   vec2eigen(pc.mid, mV);
//   vec2eigen(pc.F, mF);
//   igl::upsample(mV, mF, subV, subF, 2);

//   GEO::Mesh geoM;
//   auto [pcV, pcT] = tetmesh_from_prismcage(pc.base, pc.top, pc.F);
//   to_geogram_mesh(pcV, pcT, geoM);
//   auto tetaabb = GEO::MeshCellsAABB(geoM, false);

//   igl::embree::EmbreeIntersector embree;
//   embree.init(Eigen::MatrixXf(pc.ref.V.template cast<float>()), Eigen::MatrixXi(pc.ref.F), true);

//   RowMatd result(subV.rows(), 3);  // fid, u,v
//   for (int i = 0; i < subV.rows(); i++) {
//     spdlog::trace("subV.row({}) {}", i, subV.row(i));
//     GEO::vec3 p(subV(i, 0), subV(i, 1), subV(i, 2));
//     auto tetid = tetaabb.containing_tet(p);
//     assert(tetid != GEO::MeshCellsAABB::NO_TET);
//     auto fid = tetid / 3;
//     auto pillar_id = tetid % 3;
//     spdlog::trace("pr {} pi {}", fid, pillar_id);

//     auto [v0, v1, v2] = pc.F[fid];

//     Vec3d abt(0, 0, 0);  // xyz in the standard prism
//     {                    // projection to standard abt
//       Eigen::RowVector4d bary;
//       auto& pctet = pcT[tetid];
//       igl::barycentric_coordinates(subV.row(i), pcV[pctet[0]], pcV[pctet[1]],
//                                    pcV[pctet[2]], pcV[pctet[3]], bary);
//       spdlog::trace("Tet Bary {}", bary);
//       auto tetc = (v1 > v2 ? TETRA_SPLIT_A : TETRA_SPLIT_B)[pillar_id];
//       for (int j = 0; j < 4; j++) abt += bary[j] * STANDARD_PRISM[tetc[j]];
//     }
//     Vec3d facet_bary(1 - abt[0] - abt[1], abt[0], abt[1]);
//     auto inter = fiber_hit_ref(pc.base, pc.top, embree, v0, v1, v2, facet_bary);

//     result.row(i) = inter;
//   }
// }