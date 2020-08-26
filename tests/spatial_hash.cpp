#include <doctest.h>
#include <geogram/basic/geometry.h>
#include <igl/avg_edge_length.h>
#include <igl/read_triangle_mesh.h>
#include <prism/spatial-hash/AABB_hash.hpp>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/spdlog.h>
#include <numeric>

TEST_CASE("spatial hash") {
  RowMatd V;
  RowMati F;
  igl::read_triangle_mesh("../tests/data/bunny.off", V, F);
  put_in_unit_box(V);
  double avg_len = igl::avg_edge_length(V, F);
  Vec3d max_bnd = V.colwise().maxCoeff();
  Vec3d min_bnd = V.colwise().minCoeff();
  std::vector<Vec3d> vecV;
  std::vector<Vec3i> vecF;
  eigen2vec(V, vecV);
  eigen2vec(F, vecF);
  prism::HashGrid hg(min_bnd, max_bnd, avg_len);
  std::vector<int> frange(vecF.size());
  std::iota(frange.begin(), frange.end(), 0);
  hg.insert_triangles(vecV, vecF, frange);
  Eigen::Matrix3d local;
  for (auto k : {0, 1, 2})
    local.row(k) = V.row(F(0, k));
  auto aabb_min = local.colwise().minCoeff();
  auto aabb_max = local.colwise().maxCoeff();
  std::set<int> q;
  hg.query(aabb_min, aabb_max, q);
  CHECK(q.size() == 21);
}