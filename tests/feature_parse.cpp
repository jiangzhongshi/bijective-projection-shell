#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/local_operations/remesh_pass.hpp"
#include "test_common.hpp"
#include <doctest.h>
#include <fstream>
#include <igl/read_triangle_mesh.h>
#include <prism/feature_utils.hpp>
#include <prism/geogram/geogram_utils.hpp>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/spdlog.h>
#include <string>


TEST_CASE("feature datagen") {
  prism::geo::init_geogram();
  std::string filepath("/home/zhongshi/data/gao_hex/cad/block_input_tri");
  if (true) {
    RowMatd V;
    RowMati F;
    igl::read_triangle_mesh(filepath + ".obj", V, F);
    put_in_unit_box(V);
    Eigen::VectorXi feature_corners;
    RowMati feature_edges;
    prism::read_feature_graph(filepath + ".fgraph", feature_corners,
                              feature_edges);
    PrismCage pc(V, F, -1 /*disable bevel*/, /*initial_thick=*/0.01,
                 PrismCage::SeparateType::kShell,
                 feature_corners, feature_edges);
    std::vector<std::list<int>> chains;
    prism::split_feature_graph(pc.feature_corners, pc.feature_edges,
                               pc.mid.size(), chains);

    std::vector<std::set<int>> feature_region_segments;
    prism::feature_chain_region_segments(pc.ref.V, pc.ref.F, chains,
                                         feature_region_segments);
    for (int i = 0; i < chains.size(); i++) {
      auto &left = feature_region_segments[i * 2],
           right = feature_region_segments[i * 2 + 1];
      for (auto f : left) {
        std::set<int> s;
        std::set_difference(pc.track_ref[f].begin(), pc.track_ref[f].end(),
                            right.begin(), right.end(),
                            std::inserter(s, s.begin()));
        pc.track_ref[f] = std::move(s);
      }
      for (auto f : right) {
        std::set<int> s;
        std::set_difference(pc.track_ref[f].begin(), pc.track_ref[f].end(),
                            left.begin(), left.end(),
                            std::inserter(s, s.begin()));
        pc.track_ref[f] = std::move(s);
      }
    }

    pc.serialize("block_init.h5");
  }
  PrismCage pc("block_init.h5");
  prism::local::RemeshOptions options(pc.mid.size(), 0.5);
  options.distortion_bound = 1e-5;
  options.collapse_quality_threshold = 1e3;
  options.target_thickness = 0.2;
  options.target_adjustment.resize(pc.mid.size(), 1);

  options.collapse_quality_threshold = 30;
  options.parallel = false;
  //  int col = prism::local::wildcollapse_pass(pc, options);
  for (int i = 0; i < 2; i++) {
    prism::local::localsmooth_pass(pc, options);
    prism::local::wildflip_pass(pc, options);
    prism::local::wildcollapse_pass(pc, options);
  }
  pc.serialize("block_collapse.h5");
}

bool cage_is_away_from_ref(PrismCage &pc) {
  pc.ref.aabb.reset(new prism::geogram::AABB(pc.ref.V, pc.ref.F));
  for (auto [v0, v1, v2] : pc.F) {
    // singular edge does not exist. Bevel always split it aggressively.
    assert(!(v1 < pc.ref.aabb->num_freeze && v2 < pc.ref.aabb->num_freeze));
    if (pc.ref.aabb->intersects_triangle(
            {pc.base[v0], pc.base[v1], pc.base[v2]},
            v0 < pc.ref.aabb->num_freeze)) {
      spdlog::dump_backtrace();
      spdlog::error("Intersect Base at {} {} {}", v0, v1, v2);
      exit(1);
      return false;
    }
    if (pc.ref.aabb->intersects_triangle({pc.top[v0], pc.top[v1], pc.top[v2]},
                                         v0 < pc.ref.aabb->num_freeze)) {
      spdlog::error("Intersect Top at {} {} {}", v0, v1, v2);
      exit(1);
      return false;
    }
  }
  return true;
}

namespace prism::local {
int featurecollapse_pass(PrismCage &, RemeshOptions &,
                         std::vector<std::set<int>> &);
}

TEST_CASE("feature simplify") {
  prism::geo::init_geogram();
  PrismCage pc("../build_profile/block_collapse.h5");
  prism::local::RemeshOptions options(pc.mid.size(), 0.5);

  options.distortion_bound = 1e-5;
  options.target_thickness = 0.2;
  options.target_adjustment.resize(pc.mid.size(), 1);
  options.collapse_quality_threshold = 30;
  options.parallel=true;

  std::vector<std::list<int>> chains;
  prism::split_feature_graph(pc.feature_corners, pc.feature_edges,
                             pc.mid.size(), chains);

  std::vector<std::set<int>> feature_region_segments;
  prism::feature_chain_region_segments(pc.ref.V, pc.ref.F, chains,
                                       feature_region_segments);
  std::map<std::pair<int, int>, std::pair<int, std::vector<int>>> meta_edges;
  for (auto cid = 0; cid < chains.size(); cid++) {
    auto &c = chains[cid];
    for (auto it = std::next(c.begin()); it != c.end(); it++) {
      auto pit = std::prev(it);
      meta_edges[{*pit, *it}] = {cid, {*pit, *it}};
    }
  }

  pc.meta_edges = std::move(meta_edges);
  spdlog::set_level(spdlog::level::info);
  prism::local::featurecollapse_pass(pc, options, feature_region_segments);

  prism::local::featurecollapse_pass(pc, options, feature_region_segments);
  for (int i = 0; i < 20; i++) {
    prism::local::localsmooth_pass(pc, options);
    prism::local::wildflip_pass(pc, options);
    prism::local::wildcollapse_pass(pc, options);
    prism::local::featurecollapse_pass(pc, options, feature_region_segments);
    cage_is_away_from_ref(pc);
  }
  spdlog::info("meta_edges\n {}", pc.meta_edges);
  pc.serialize("block_feature.h5");
}
