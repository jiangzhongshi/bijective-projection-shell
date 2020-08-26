#include "test_common.hpp"
#include <doctest.h>
#include <prism/local_operations/local_mesh_edit.hpp>
#include <spdlog/spdlog.h>

#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/local_operations/remesh_pass.hpp"
#include "prism/predicates/positive_prism_volume_12.hpp"
#include "prism/spatial-hash/AABB_hash.hpp"
#include <igl/read_triangle_mesh.h>
#include "prism/geogram/geogram_utils.hpp"

void zoom_out(PrismCage &pc, double eps) {
  for (int i = 0; i < pc.mid.size(); i++) {
    pc.base[i] = (1 - eps) * pc.mid[i] + pc.base[i] * eps;
    pc.top[i] = (1 - eps) * pc.mid[i] + pc.top[i] * eps;
  }
}

TEST_CASE("bunny collapse") {
  prism::geo::init_geogram();
  spdlog::set_level(spdlog::level::info);
  if (false) {
    RowMatd V;
    RowMati F;
    igl::read_triangle_mesh(
        "/home/zhongshi/Workspace/libigl/tutorial/data/bunny.off", V, F);
    put_in_unit_box(V);
    PrismCage pc(V, F, 0.2, 0.1);
    pc.serialize("bunny_temp.h5");
  }
  PrismCage pc("../build_debug/bunny_temp.h5");
  // pc.base_grid = std::make_shared<prism::HashGrid>(pc.base, pc.F);
  // pc.top_grid = std::make_shared<prism:: HashGrid>(pc.top, pc.F);
  prism::local::RemeshOptions options(pc.mid.size(), 0.5);
  options.distortion_bound = 1e-5;
  options.target_adjustment.resize(pc.mid.size(), 1);

  options.collapse_quality_threshold = 30;
  options.parallel = false;

  int col = prism::local::wildcollapse_pass(pc, options);
  // CHECK(col == 6221);
}
