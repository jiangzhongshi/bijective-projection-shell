#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/read_triangle_mesh.h>
#include <igl/upsample.h>
#include <igl/write_triangle_mesh.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <highfive/H5Easy.hpp>
#include <prism/cgal/polyhedron_self_intersect.hpp>
#include <prism/energy/prism_quality.hpp>
#include <prism/local_operations/section_remesh.hpp>

#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
extern "C" {
size_t getPeakRSS();
}

double total_energy(const std::vector<Vec3d>& V, const std::vector<Vec3i>& F);

void remesh_in_shell(const PrismCage& pc, std::vector<Vec3d>& V,
                     std::vector<Vec3i>& F, double targ_len, RowMatd& mV,
                     RowMati& mF) {
  RowMatd mT, mB;
  vec2eigen(pc.base, mB);
  vec2eigen(pc.top, mT);
  vec2eigen(pc.F, mF);

  std::vector<std::set<int>> track_to_prism;
  if (V.size() == 0) {
    eigen2vec(pc.ref.V, V);
    eigen2vec(pc.ref.F, F);
    track_to_prism.resize(F.size());
    for (int p = 0; p < pc.track_ref.size(); p++) {
    for (auto t : pc.track_ref[p]) track_to_prism[t].insert(p);
  }
  }
  else {
    spdlog::error("need to reinitialize track");
    exit(1);
  }

  spdlog::enable_backtrace(100);


  auto tree_B = prism::geogram::AABB(mB, mF);
  auto tree_T = prism::geogram::AABB(mT, mF);
  std::vector<double> target_adjustment(V.size(), 1);

  prism::section::RemeshOptions option(V.size(), targ_len);
  option.distortion_bound = 1e-5;
  option.collapse_improve_quality = false;
  total_energy(V, F);
  for (int i = 0; i < 500; i++) {
    int col = prism::section::wildcollapse_pass(
        pc, tree_B, tree_T, option, V, F, track_to_prism, target_adjustment);
    for (int j = 0; j < 4; j++) {
      prism::section::wildflip_pass(pc, tree_B, tree_T, option, V, F,
                                    track_to_prism);
      prism::section::localsmooth_pass(pc, tree_B, tree_T, option, V, F,
                                       track_to_prism);
    }
    total_energy(V, F);
    if (col < 1) break;
  }
  spdlog::info("Finish first collapse pass. Proceed to quality improvement");
  option.split_improve_quality = false;
  option.collapse_improve_quality = true;
  for (int i = 0; i < 500; i++) {
    spdlog::info("{} Start Split, #V {}", i, V.size());
    prism::section::wildsplit_pass(pc, tree_B, tree_T, option, V, F,
                                   track_to_prism, target_adjustment);
    total_energy(V, F);
    for (int j = 0; j < 3; j++) {
      prism::section::wildflip_pass(pc, tree_B, tree_T, option, V, F,
                                    track_to_prism);
      prism::section::localsmooth_pass(pc, tree_B, tree_T, option, V, F,
                                       track_to_prism);
    }
    total_energy(V, F);
    spdlog::info("{} Start Collapse, #V {}", i, V.size());

    int col = prism::section::wildcollapse_pass(
        pc, tree_B, tree_T, option, V, F, track_to_prism, target_adjustment);
    total_energy(V, F);
    for (int j = 0; j < 3; j++) {
      total_energy(V, F);
      prism::section::wildflip_pass(pc, tree_B, tree_T, option, V, F,
                                    track_to_prism);
      total_energy(V, F);
      prism::section::localsmooth_pass(pc, tree_B, tree_T, option, V, F,
                                       track_to_prism);
    }
    prism::section::wildcollapse_pass(
        pc, tree_B, tree_T, option, V, F, track_to_prism, target_adjustment);
    total_energy(V, F);
    if (col < 1) break;
    vec2eigen(V, mV);
    vec2eigen(F, mF);
    igl::write_triangle_mesh("temp.ply", mV, mF);
  }

  vec2eigen(V, mV);
  vec2eigen(F, mF);
}

void section_pipeline(std::string filename, std::string sec_input,
                      std::string ser_file, double target_edge = 0.1) {
  PrismCage pc(filename);
  RowMatd mV;
  RowMati mF;
  std::vector<Vec3d> Vin;
  std::vector<Vec3i> Fin;
  if (!sec_input.empty()) {
    spdlog::info("loading sec {}", sec_input);
    igl::read_triangle_mesh(sec_input, mV, mF);
    put_in_unit_box(mV);
    eigen2vec(mV, Vin);
    eigen2vec(mF, Fin);
  }
  remesh_in_shell(pc, Vin, Fin, target_edge, mV, mF);
  spdlog::info("Peak Memory {}", getPeakRSS());
  igl::writePLY(ser_file, mV, mF);  // default binary
}