#include <doctest.h>

#include <geogram/mesh/mesh_AABB.h>
#include <igl/exact_geodesic.h>
#include <igl/heat_geodesics.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>
#include <prism/cage_utils.hpp>
#include <prism/common.hpp>
#include <prism/geogram/geogram_utils.hpp>
#include <prism/phong/projection.hpp>
#include <highfive/H5Easy.hpp>
#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
#include "test_common.hpp"


TEST_CASE("Heat Geodesic") {
  RowMatd V;
  RowMati F;
  {
    Eigen::VectorXi SVI, SVJ;
    igl::read_triangle_mesh(
        "/home/zhongshi/data/Thingi10K/raw_meshes/121868.stl", V, F);
    Eigen::MatrixXd temp_V = V;  // for STL file
    igl::remove_duplicate_vertices(temp_V, 0, V, SVI, SVJ);
    for (int i = 0; i < F.rows(); i++)
      for (int j : {0, 1, 2}) F(i, j) = SVJ[F(i, j)];

    spdlog::info("V={}, F={}", V.rows(), F.rows());
    put_in_unit_box(V);
  }

  igl::HeatGeodesicsData<double> data;
  double t = std::pow(igl::avg_edge_length(V, F), 2);
  if (!igl::heat_geodesics_precompute(V, F, t, data)) {
    std::cerr << "Error: heat_geodesics_precompute failed." << std::endl;
    exit(EXIT_FAILURE);
  };
  Eigen::VectorXd D = Eigen::VectorXd::Zero(data.Grad.cols());
  D(114) = 1;
  igl::heat_geodesics_solve(data, (Eigen::VectorXi(1, 1) << 114).finished(), D);
  Eigen::VectorXi VT;
  VT.setLinSpaced(V.rows(), 0, V.rows() - 1);
  Eigen::VectorXd DE;
  igl::exact_geodesic(V, F, (Eigen::VectorXi(1, 1) << 114).finished(),
                      Eigen::VectorXi(), VT, Eigen::VectorXi(), DE);
  H5Easy::File file("../tests/data/heat.h5", H5Easy::File::Overwrite);
  H5Easy::dump(file, "DH", D);
  H5Easy::dump(file, "DE", DE);
  H5Easy::dump(file, "V", V);
  H5Easy::dump(file, "F", F);
}

void heat_geodesic(const RowMatd& V, const RowMati& F, Vec3d& target,
                   Eigen::VectorXd& D) {
  igl::HeatGeodesicsData<double> data;
  double t = std::pow(igl::avg_edge_length(V, F), 2);
  igl::heat_geodesics_precompute(V, F, t, data);
  int vid = -1;
  double min_norm = 1;
  for (int i = 0; i < V.rows(); i++) {
    if ((V.row(i) - target).norm() < min_norm) {
      min_norm = (V.row(i) - target).norm();
      vid = i;
    }
  }
  D = Eigen::VectorXd::Zero(data.Grad.cols());
  D(vid) = 1;
  igl::heat_geodesics_solve(data, (Eigen::VectorXi(1, 1) << vid).finished(), D);
}
void exact_geodesic(RowMatd& V, RowMatd& F, Vec3d& target, Eigen::VectorXd& D) {
  int vid = -1;
  double min_norm = 1;
  for (int i = 0; i < V.rows(); i++) {
    if ((V.row(i) - target).norm() < min_norm) {
      min_norm = (V.row(i) - target).norm();
      vid = i;
    }
  }

  Eigen::VectorXi VT;
  VT.setLinSpaced(V.rows(), 0, V.rows() - 1);
  igl::exact_geodesic(V, F, (Eigen::VectorXi(1, 1) << vid).finished(),
                      Eigen::VectorXi(), VT, Eigen::VectorXi(), D);
}

#include <igl/Hit.h>
bool project_to_proxy_mesh(const PrismCage& pc, const prism::geogram::AABB& pxtree,
                           int prism_id, const Vec3d& spatial, igl::Hit& hit) {
  std::array<double, 3> tuple;
  auto [v0, v1, v2] = pc.F[prism_id];
  bool bottom = false;
  std::array<Vec3d, 6> stackV = {pc.base[v0], pc.base[v1], pc.base[v2],
                                 pc.mid[v0],  pc.mid[v1],  pc.mid[v2]};
  bool hitbottom =
      prism::phong::phong_projection(stackV, spatial, v1 > v2, tuple);
  stackV = {pc.mid[v0], pc.mid[v1], pc.mid[v2],
            pc.top[v0], pc.top[v1], pc.top[v2]};
  bool hittop = prism::phong::phong_projection(stackV, spatial, v1 > v2, tuple);
  if (!hitbottom && !hittop) {
    spdlog::error("No proj");
    return false;
  }
  std::array<Vec3d, 4> endpoints0;
  prism::phong::fiber_endpoints({pc.base[v0], pc.base[v1], pc.base[v2],
                                 pc.mid[v0], pc.mid[v1], pc.mid[v2]},
                                v1 > v2, tuple[0], tuple[1], endpoints0);
  for (int i = 0; i < 3; i++)
    if (pxtree.segment_hit(endpoints0[i], endpoints0[i + 1], hit)) return true;
  std::array<Vec3d, 4> endpoints1;
  prism::phong::fiber_endpoints(
      {pc.mid[v0], pc.mid[v1], pc.mid[v2], pc.top[v0], pc.top[v1], pc.top[v2]},
      v1 > v2, tuple[0], tuple[1], endpoints1);
  for (int i = 0; i < 3; i++)
    if (pxtree.segment_hit(endpoints1[i], endpoints1[i + 1], hit)) return true;
  return false;
}
/*python script for rendering
with h5py.File('../buildr/heat.h5') as fp:
    V0,F0,D0,E0 = map(lambda x:fp[x][()], ['V','F','DH','DE'])
with h5py.File('../buildr/heat2.h5') as fp:
    V1,F1,nutD,transfered = map(lambda x:fp[x][()], ['V','F','DH','qD'])

plt = None
V0 = V0[:,[1,2,0]]
V1 = V1[:,[1,2,0]]
sh = dict(width=1000,height=1000)
vw = mp.Viewer(sh)
vw.add_mesh(V0,F0,c=D0,shading=sh)
vw.add_edges(*igl.isolines(V0,F0,D0,10))
plt = mp.Subplot(plt, vw,s=[2,2,0])

vw = mp.Viewer(sh)
vw.add_mesh(V0,F0,c=E0,shading=sh)
vw.add_edges(*igl.isolines(V0,F0,E0,10))
plt = mp.Subplot(plt, vw,s=[2,2,1])

vw = mp.Viewer(sh)
vw.add_mesh(V1,F1,c=nutD,shading=sh)
vw.add_edges(*igl.isolines(V1,F1,nutD,10))
plt = mp.Subplot(plt, vw,s=[2,2,2])

vw = mp.Viewer(sh)
vw.add_mesh(V0,F0,c=transfered,shading=sh)
vw.add_edges(*igl.isolines(V0,F0,transfered,10))
plt = mp.Subplot(plt, vw,s=[2,2,3])
*/
TEST_CASE("NutHeat") {
  RowMatd pxV;
  RowMati pxF;
  igl::read_triangle_mesh("../tests/data/heat_swirl.obj", pxV, pxF);
  PrismCage pc("../tests/data/121868.stl.h5");
  Eigen::VectorXd pxD, DE;
  Vec3d target(-0.08108107, -0.5781373, -0.74871287);
  heat_geodesic(pxV, pxF, target, pxD);
  // pxD.resize(pxV.rows());
  // heat_geodesic(pxV, pxF, target, DE);

  // transfer: for each point on the raw mesh.
  H5Easy::File file0("../tests/data/heat.h5", H5Easy::File::ReadOnly);
  RowMatd meshV = H5Easy::load<RowMatd>(file0, "V");
  RowMati meshF = H5Easy::load<RowMati>(file0, "F");

  std::vector<Vec3d> tetV;
  std::vector<Vec4i> tetT;
  prism::cage_utils::tetmesh_from_prismcage(pc.base, pc.mid, pc.top, pc.F, 0, tetV,
                                            tetT);
  RowMatd mtetV;
  RowMati mtetT;
  vec2eigen(tetV, mtetV);
  vec2eigen(tetT, mtetT);
  prism::geo::init_geogram();
  GEO::Mesh geo_tet;
  prism::geo::to_geogram_mesh(mtetV, mtetT, geo_tet);
  GEO::MeshCellsAABB tetaabb(geo_tet, false);
  prism::geogram::AABB pxtree(pxV, pxF);

  Eigen::VectorXd queryD(meshV.rows());
  for (int i = 0; i < meshV.rows(); i++) {
    Vec3d spatial = meshV.row(i);
    auto tet_id =
        tetaabb.containing_tet(GEO::vec3(spatial[0], spatial[1], spatial[2]));
    auto prism_id = tet_id / 6;
    igl::Hit hit{-1, -1, -1, -1, -1};
    if (!project_to_proxy_mesh(pc, pxtree, prism_id, spatial, hit)) {
      spdlog::error("Failed {}", i);
      exit(1);
    };
    queryD[i] = pxD[pxF(hit.id, 0)] * (1 - hit.u - hit.v) +
                pxD[pxF(hit.id, 1)] * (hit.u) + pxD[pxF(hit.id, 2)] * (hit.v);
  }

  H5Easy::File file("../tests/data/heat2.h5", H5Easy::File::Overwrite);
  H5Easy::dump(file, "DH", pxD);
  H5Easy::dump(file, "qD", queryD);
  H5Easy::dump(file, "V", pxV);
  H5Easy::dump(file, "F", pxF);
}

TEST_CASE("Wheel Elasticity") {
  RowMatd pxV;
  RowMati pxF;
  H5Easy::File file("../buildr/wheel_stretch.h5", H5Easy::File::ReadOnly);
  RowMatd pxD;
  RowMatd meshV;
  RowMati meshF;
  pxV = H5Easy::load<RowMatd>(file, "tV");
  pxD = H5Easy::load<RowMatd>(file, "tD");
  meshV = H5Easy::load<RowMatd>(file, "mV");
  pxF = H5Easy::load<RowMati>(file, "tF");
  meshF = H5Easy::load<RowMati>(file, "mF");
  PrismCage pc("../tests/data/1517923.obj.h5");

  std::vector<Vec3d> tetV;
  std::vector<Vec4i> tetT;
  prism::cage_utils::tetmesh_from_prismcage(pc.base, pc.mid, pc.top, pc.F, 0, tetV,
                                            tetT);
  RowMatd mtetV;
  RowMati mtetT;
  vec2eigen(tetV, mtetV);
  vec2eigen(tetT, mtetT);
  prism::geo::init_geogram();
  GEO::Mesh geo_tet;
  prism::geo::to_geogram_mesh(mtetV, mtetT, geo_tet);
  GEO::MeshCellsAABB tetaabb(geo_tet, false);
  prism::geogram::AABB pxtree(pxV, pxF);
  RowMatd queryD(meshV.rows(),3);
  for (int i = 0; i < meshV.rows(); i++) {
    Vec3d spatial = meshV.row(i);
    auto tet_id =
        tetaabb.containing_tet(GEO::vec3(spatial[0], spatial[1], spatial[2]));
    auto prism_id = tet_id / 6;
    igl::Hit hit{-1, -1, -1, -1, -1};
    if (!project_to_proxy_mesh(pc, pxtree, prism_id, spatial, hit)) {
      spdlog::error("Failed {}", i);
      exit(1);
    };
    queryD.row(i) = pxD.row(pxF(hit.id, 0)) * (1 - hit.u - hit.v) +
                pxD.row(pxF(hit.id, 1)) * (hit.u) + pxD.row(pxF(hit.id, 2)) * (hit.v);
  }
  igl::write_triangle_mesh("../tests/data/wheel_transfered.obj",meshV+queryD, meshF);
  // H5Easy::File file("../tests/data/heat2.h5", H5Easy::File::Overwrite);
  // H5Easy::dump(file, "DH", pxD);
  // H5Easy::dump(file, "qD", queryD);
  // H5Easy::dump(file, "V", pxV);
  // H5Easy::dump(file, "F", pxF);
}