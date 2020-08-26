#include <igl/per_vertex_normals.h>
#include <igl/ray_mesh_intersect.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>
#include <prism/common.hpp>
#include <prism/energy/map_distortion.hpp>
#include <prism/local_operations/section_remesh.hpp>
#include <prism/phong/projection.hpp>
#include <doctest.h>
#include <highfive/H5Easy.hpp>
#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/local_operations/remesh_pass.hpp"

constexpr auto polyline_hit = [](auto& pc, auto prism_id, auto u, auto v,
                                 auto& hit, auto& dir) {
  auto [v0, v1, v2] = pc.F[prism_id];
  std::array<Vec3d, 4> endpoints0;
  prism::phong::fiber_endpoints({pc.base[v0], pc.base[v1], pc.base[v2],
                                 pc.mid[v0], pc.mid[v1], pc.mid[v2]},
                                v1 > v2, u, v, endpoints0);
  for (int i = 0; i < 3; i++) {
    auto inter =
        pc.ref.aabb->segment_hit(endpoints0[i], endpoints0[i + 1], hit);
    if (inter) {
      dir = endpoints0[i + 1] - endpoints0[i];
      return true;
    }
  }
  std::array<Vec3d, 4> endpoints1;
  prism::phong::fiber_endpoints(
      {pc.mid[v0], pc.mid[v1], pc.mid[v2], pc.top[v0], pc.top[v1], pc.top[v2]},
      v1 > v2, u, v, endpoints1);
  for (int i = 0; i < 3; i++) {
    auto inter =
        pc.ref.aabb->segment_hit(endpoints1[i], endpoints1[i + 1], hit);
    if (inter) {
      dir = endpoints1[i + 1] - endpoints1[i];
      return true;
    }
  }
  return false;
};
#include <igl/upsample.h>

TEST_CASE("curve fit") {
  PrismCage pc("../buildr/armadillo.obj.h5");

  spdlog::info("Begin");
  std::vector<Eigen::Matrix<double, 1, 6>> uvt;
  std::vector<Vec3d> ambient_normal;
  std::vector<Eigen::Matrix<double, 1, 9>> frames;
  std::vector<Vec3d> all_hit_points;
  std::vector<int> prism_ptr{0};
  RowMatd refVN;
  igl::per_vertex_normals(pc.ref.V, pc.ref.F, refVN);
  RowMatd unitV(3, 2), usV;
  unitV << 0, 0, 1, 0, 0, 1;
  RowMati unitF(1, 3), usF;
  unitF << 0, 1, 2;
  igl::upsample(unitV, unitF, usV, usF, 5);
  for (int prism_id = 0; prism_id < pc.F.size(); prism_id++) {
    auto trackee = pc.track_ref[prism_id];
 
    std::array<Vec3d, 3> B, N, M;
    for (auto i : {0, 1, 2}) {
      M[i] = pc.mid[pc.F[prism_id][i]];
      B[i] = pc.base[pc.F[prism_id][i]];
      N[i] = pc.top[pc.F[prism_id][i]] - B[i];
    }
    for (int i = 0; i < usV.rows(); i++) {
      auto u = usV(i, 0), v = usV(i, 1);
      igl::Hit hit{-1, -1, -1, -1, 2};
      Vec3d dir;
      if (!polyline_hit(pc, prism_id, u, v, hit, dir)) {
        spdlog::warn("failed hit");
        continue;
      }
      if (hit.t > 1) continue;
      Vec3d hitpoint, hitnormal;
      auto ru = hit.u, rv = hit.v;
      auto face = pc.ref.F.row(hit.id);
      hitnormal = refVN.row(face[0]) * (1 - ru - rv) + refVN.row(face[1]) * ru +
                  refVN.row(face[2]) * (rv);
      hitnormal.normalize();
      hitpoint = pc.ref.V.row(face[0]) * (1 - ru - rv) +
                 pc.ref.V.row(face[1]) * ru + pc.ref.V.row(face[2]) * (rv);
      Eigen::Matrix<double, 1, 6> uvtn;
      uvtn << u, v, hit.t, hitpoint[0], hitpoint[1], hitpoint[2];
      uvt.emplace_back(uvtn);
    }
    prism_ptr.push_back(uvt.size());
  }
  spdlog::info("End, size {}, uvt {}", pc.F.size(), uvt.size());
  H5Easy::File file("armadillo_uvt.h5", H5Easy::File::Overwrite);
  RowMatd muvt = Eigen::Map<RowMatd>(uvt[0].data(), uvt.size(), 6);
  RowMatd mhitnormals =
      Eigen::Map<RowMatd>(ambient_normal[0].data(), ambient_normal.size(), 3);
  RowMatd mhitpoints =
      Eigen::Map<RowMatd>(all_hit_points[0].data(), all_hit_points.size(), 3);
  H5Easy::dump(file, "uvt", muvt);
  H5Easy::dump(file, "prism_ptr", prism_ptr);
}