#include <prism/common.hpp>
#include <doctest.h>

#include <igl/read_triangle_mesh.h>
#include <igl/volume.h>
#include <spdlog/spdlog.h>
#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/local_operations/remesh_pass.hpp"

#include "test_common.hpp"

#include <prism/energy/prism_quality.hpp>
constexpr auto total_energy = [](auto& V, auto& F) {
  std::set<int> low_quality_vertices;
  double total_quality = 0;
  double max_quality = 0;
  for (auto [v0, v1, v2] : F) {
    auto q = prism::energy::triangle_quality({V[v0], V[v1], V[v2]});
    total_quality += q;
    max_quality = std::max(max_quality, q);
  }

  spdlog::info("Total Q {} number {}, divide {}, max {}", total_quality, F.size(),
               total_quality / F.size(), max_quality);
  if (max_quality < 10){spdlog::info("SUCCESS"); exit(0);}
  return max_quality;
};


TEST_CASE("cubes doo sabin") {
  spdlog::set_level(spdlog::level::info);
  RowMatd V;
  RowMati F;
  igl::read_triangle_mesh("/home/zhongshi/Workspace/libigl/tutorial/data/camelhead.off", V, F);
  PrismCage pc(V, F, 0.2);
  pc.target_adjustment.resize(pc.mid.size(), 1);
  prism::local::RemeshOptions options(pc.mid.size(), 0.1);
  options.distortion_bound = 1e2;
  for (int i = 0; i < 5; i++) {
    if (i % 2 == 0) options.collapse_improve_quality = true;
    spdlog::info("V {} F{}", pc.mid.size(), pc.F.size());
    prism::local::localsmooth_pass(pc, options);
    pc.serialize(std::to_string(i) + "It_cube_0smooth_ds.h5");
    prism::local::wildcollapse_pass(pc, options);
    pc.serialize(std::to_string(i) + "It_cube_1collapse_ds.h5");
    prism::local::wildsplit_pass(pc, options);
    pc.serialize(std::to_string(i) + "It_cube_2split_ds.h5");
    prism::local::wildflip_pass(pc, options);
    prism::local::localsmooth_pass(pc, options);
    prism::local::wildflip_pass(pc, options);
    pc.serialize(std::to_string(i) + "It_cube_3smooth_ds.h5");
    spdlog::info("Outer iteration {} Finished", i);
  }
}


TEST_CASE("manual wheel") {
  spdlog::set_level(spdlog::level::info);
  RowMatd V;
  RowMati F;
  igl::read_triangle_mesh("../tests/data/1517923.obj", V, F);
  PrismCage pc(V, F, 0.2);
  pc.target_adjustment.resize(pc.mid.size(), 1);
  prism::local::RemeshOptions options(pc.mid.size(), 0.1);
  options.distortion_bound = 1e5;
  pc.serialize("wheelinit.h5");
  options.collapse_improve_quality = true;
  for (int i = 0; i < 50; i++) {
    spdlog::info("V {} F{}", pc.mid.size(), pc.F.size());
    prism::local::localsmooth_pass(pc, options);
    pc.serialize(std::to_string(i) + "It_wheel_0smooth_ds.h5");
    total_energy(pc.mid, pc.F);
    prism::local::wildcollapse_pass(pc, options);
    pc.serialize(std::to_string(i) + "It_wheel_1collapse_ds.h5");
    total_energy(pc.mid, pc.F);
    prism::local::wildflip_pass(pc, options);
    prism::local::localsmooth_pass(pc, options);
    pc.serialize(std::to_string(i) + "It_wheel_3smooth_ds.h5");
    total_energy(pc.mid, pc.F);
    spdlog::info("Outer iteration {} Finished", i);
  }
}


TEST_CASE("refine wheel") {
  RowMatd V;
  RowMati F;
  PrismCage pc("49It_wheel_3smooth_ds.h5");
  pc.target_adjustment.resize(pc.mid.size(), 1);
  prism::local::RemeshOptions options(pc.mid.size(), 0.1);
  options.distortion_bound = 1e10;
  for (int i = 0; i < 10; i++) {
    options.collapse_improve_quality = true;
    prism::local::wildcollapse_pass(pc, options);
    total_energy(pc.mid, pc.F);
    prism::local::wildflip_pass(pc, options);
    prism::local::localsmooth_pass(pc, options);
    total_energy(pc.mid, pc.F);
    pc.serialize(std::to_string(i) + "distorted_collapse.h5");
    }
}

TEST_CASE("manual vh") {
  spdlog::set_level(spdlog::level::info);
  RowMatd V;
  RowMati F;
  igl::read_triangle_mesh(
    "/home/zhongshi/data/mpz14_data/inputmodels/vh_skin.obj", V, F);
  put_in_unit_box(V);
  PrismCage pc(V, F, 0.2);
  pc.target_adjustment.resize(pc.mid.size(), 1);
  prism::local::RemeshOptions options(pc.mid.size(), 0.1);
  options.distortion_bound = 1e10;
  for (int i = 0; i < 50; i++) {
    options.collapse_improve_quality = false;
    if (i % 2 == 0) options.collapse_improve_quality = true;
    spdlog::info("V {} F{}", pc.mid.size(), pc.F.size());
    prism::local::localsmooth_pass(pc, options);
    total_energy(pc.mid, pc.F);
    prism::local::wildcollapse_pass(pc, options);
    total_energy(pc.mid, pc.F);
    prism::local::wildflip_pass(pc, options);
    prism::local::localsmooth_pass(pc, options);
    pc.serialize(std::to_string(i) + "It_vh_ds.h5");
    total_energy(pc.mid, pc.F);
    spdlog::info("Outer iteration {} Finished", i);
  }
}