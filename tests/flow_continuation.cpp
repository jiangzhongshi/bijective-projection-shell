#include <spdlog/spdlog.h>
#include <doctest.h>

#include <igl/triangle_triangle_adjacency.h>
#include <igl/volume.h>
#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/local_operations/remesh_pass.hpp"
#include "prism/predicates/positive_prism_volume_12.hpp"
#include "test_common.hpp"
#include <prism/extraction.hpp>


void reset_ref_to_current(PrismCage& pc) {
  vec2eigen(pc.mid, pc.ref.V);
  vec2eigen(pc.F, pc.ref.F);
  pc.ref.aabb = std::make_unique<prism::geogram::AABB>(pc.ref.V, pc.ref.F);
  pc.init_track();
}

TEST_CASE("blow up") {
  PrismCage pc("../buildr/armadillo.obj.h5");
  RowMatd V;
  RowMati F;
  prism::local::RemeshOptions options;
  spdlog::set_level(spdlog::level::info);
  options.collapse_improve_quality = false;
  for (int blow = 0; blow < 3; blow++) {
    vec2eigen(pc.top, V);
    vec2eigen(pc.F, F);
    pc = PrismCage(V, F, 0.25);
    spdlog::info("Start Pass #V {}", pc.mid.size());
    for (int i = 0; i < 10; i++) {
      options.sizing_field = [&blow](const Vec3d&) -> double { return 1e-1; };
      int collapsed = prism::local::wildcollapse_pass(pc, options);
      

      spdlog::info("End Collapse, #V {}", pc.mid.size());
      for (int j = 0; j < 5; j++) {
        prism::local::wildflip_pass(pc, options);
        spdlog::info("End Flip");
        prism::local::localsmooth_pass(pc, options);
        spdlog::info("End Pass");
      }

      if (collapsed < 1) break;
    }
    prism::shell_extraction(pc, false);
    std::string filename = std::to_string(blow) + "blow.h5";
    pc.serialize(filename);
  }
}
