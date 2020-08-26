
#include <igl/boundary_loop.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/writeOBJ.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <highfive/H5Easy.hpp>
#include <prism/cgal/polyhedron_self_intersect.hpp>
#include <prism/energy/map_distortion.hpp>
#include <prism/energy/prism_quality.hpp>
#include <prism/extraction.hpp>
#include <prism/geogram/geogram_utils.hpp>
#include <prism/local_operations/remesh_pass.hpp>
#include <prism/local_operations/validity_checks.hpp>

#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/predicates/positive_prism_volume_12.hpp"

extern "C" {
size_t getPeakRSS();
}

bool cage_is_away_from_ref(const PrismCage &pc) {
  for (auto [v0, v1, v2] : pc.F) {
    // singular edge does not exist. Bevel always split it aggressively.
    assert(!(v1 < pc.ref.aabb->num_freeze && v2 < pc.ref.aabb->num_freeze));
    if (pc.ref.aabb->intersects_triangle(
            {pc.base[v0], pc.base[v1], pc.base[v2]},
            v0 < pc.ref.aabb->num_freeze)) {
      spdlog::dump_backtrace();
      spdlog::debug("Intersect Base at {} {} {}", v0, v1, v2);
      return false;
    }
    if (pc.ref.aabb->intersects_triangle({pc.top[v0], pc.top[v1], pc.top[v2]},
                                         v0 < pc.ref.aabb->num_freeze)) {
      spdlog::debug("Intersect Top at {} {} {}", v0, v1, v2);
      return false;
    }
  }
  return true;
}

#include <igl/doublearea.h>
#include <igl/per_face_normals.h>
double min_dihedral_angles(const RowMatd &V, const RowMati &F) {
  RowMatd FN;
  RowMati TT, TTi;
  double minangle = 1;
  igl::per_face_normals_stable(V, F, FN);

  igl::triangle_triangle_adjacency(F, TT, TTi);
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      if (TT(i, j) == -1)
        continue;
      double a = FN.row(TT(i, j)).dot(FN.row(i));
      minangle = std::min(minangle, a);
    }
  }
  return minangle;
}

bool preconditions(const RowMatd &V, const RowMati &F,
                   const std::string &filename,
                   const std::map<std::string, bool> &controls) {
  if (F.rows() == 0) {
    spdlog::error("Precondition: Empty Mesh");
    return false;
  }
  Eigen::VectorXi BI;
  if (!igl::is_vertex_manifold(F, BI) || !igl::is_edge_manifold(F)) {
    spdlog::error("Precondition: Input {} is not manifold", filename);
    return false;
  }
  Eigen::VectorXi bnd;
  igl::boundary_loop(F, bnd);
  if (bnd.size() > 0) {
    spdlog::info("Input {} has boundary", filename);
  }

  if (controls.at("allow-intersect"))
    return true;
  Eigen::VectorXd M;
  igl::doublearea(V, F, M);
  double minarea = M.minCoeff() / 2;
  if (minarea < 1e-16) {
    spdlog::error("Precondition: Input {} area small {}", filename, minarea);
    return false;
  }
  double minangle = min_dihedral_angles(V, F);
  spdlog::info("Minimum Angle {:.18f}", minangle);
  if (minangle < -1 + 1.5230867123072755e-06) { // smaller than 0.1 degree
    spdlog::error("Precondition: Input {} flat angle", filename);
    return false;
  }
  // numerical self intersection
  prism::geogram::AABB tree(V, F, controls.at("allow-intersect") == false);
  if (tree.numerical_self_intersection(1e-10)) {
    spdlog::error("Precondition: Input {} N-self intersects", filename);
    return false;
  }
  if (prism::cgal::polyhedron_self_intersect(V, F)) {
    spdlog::error("Precondition: Input {} self intersects", filename);
    return false;
  }
  return true;
}

double total_energy(const std::vector<Vec3d> &V, const std::vector<Vec3i> &F) {
  std::set<int> low_quality_vertices;
  double total_quality = 0;
  double max_quality = 0;
  for (auto [v0, v1, v2] : F) {
    auto q = prism::energy::triangle_quality({V[v0], V[v1], V[v2]});
    total_quality += q;
    max_quality = std::max(max_quality, q);
  }

  spdlog::info("Total Q {} fnum {}, avg {}, max {}", total_quality, F.size(),
               total_quality / F.size(), max_quality);
  return max_quality;
};

double total_vol_energy(const std::vector<Vec3d> &base,
                        const std::vector<Vec3d> &mid,
                        const std::vector<Vec3d> &top,
                        const std::vector<Vec3i> &F) {
  std::set<int> low_quality_vertices;
  double total_quality = 0;
  double max_quality = 0;
  for (auto [v0, v1, v2] : F) {
    auto q = prism::energy::prism_full_quality(
                 {base[v0], base[v1], base[v2], mid[v0], mid[v1], mid[v2]}) +
             prism::energy::prism_full_quality(
                 {mid[v0], mid[v1], mid[v2], top[v0], top[v1], top[v2]});
    total_quality += q;
    max_quality = std::max(max_quality, q);
  }

  spdlog::info("Total vQ {} number {}, divide {}, max {}", total_quality,
               F.size(), total_quality / F.size(), max_quality);
  return max_quality;
};

int remesh_schedule(PrismCage &pc, prism::local::RemeshOptions &options,
                    std::string ser_file) {
  options.split_improve_quality = true;
  double old_quality = total_energy(pc.mid, pc.F);
  double cur_quality = 20;
  for (int i = 0; i < 500; i++) {
    // prism::local::wildflip_pass(pc, options);
    // total_energy(pc.mid, pc.F);
    // prism::local::localsmooth_pass(pc, options);
    // total_energy(pc.mid, pc.F);

    spdlog::info("{} Start Collapse, #V {}", i, pc.mid.size());
    int collapse_count = prism::local::wildcollapse_pass(pc, options);
    spdlog::info("End Collapse, #V {}", pc.mid.size());
    for (int j = 0; j < 2; j++) {
      prism::local::wildflip_pass(pc, options);
      total_energy(pc.mid, pc.F);
      prism::local::localsmooth_pass(pc, options);
      total_energy(pc.mid, pc.F);
      spdlog::info("End Flip+Smooth");
    }
    cur_quality = total_energy(pc.mid, pc.F);
    if (cur_quality < 4.1) {
      spdlog::info("Success Finished");
      break;
      // return 0;
    }
    if (collapse_count <= 1e-4 * pc.mid.size()) {
      spdlog::info("Converged {}-> {}, stop", old_quality, cur_quality);
      break;
    }
    old_quality = cur_quality;
  }
  for (int i = 0; i < 20; i++) {
    prism::local::wildflip_pass(pc, options);
    prism::local::localsmooth_pass(pc, options);
    spdlog::info("End Flip+Smooth");

    cur_quality = total_energy(pc.mid, pc.F);

    if (cur_quality > 0.99 * old_quality)
      break;
    old_quality = cur_quality;
  }
  return 0;
}

void shell_pipeline(std::string filename, std::string ser_file,
                    const std::map<std::string, bool> &controls,
                    double edge_bb_ratio, double angle_distortion,
                    double target_thickness, double init_thick) {
  prism::geo::init_geogram();
  spdlog::info("Input file {}", filename);
  std::unique_ptr<PrismCage> pc;
  if (std::filesystem::path(filename).extension() == ".init" ||
      std::filesystem::path(filename).extension() == ".h5") {
    pc.reset(new PrismCage(filename));
  } else {
    RowMatd V;
    RowMati F;
    {
      Eigen::VectorXi SVI, SVJ;
      igl::read_triangle_mesh(filename, V, F);
      RowMatd temp_V = V; // for STL file
      igl::remove_duplicate_vertices(temp_V, 0, V, SVI, SVJ);
      for (int i = 0; i < F.rows(); i++)
        for (int j : {0, 1, 2})
          F(i, j) = SVJ[F(i, j)];

      spdlog::info("V={}, F={}", V.rows(), F.rows());
      put_in_unit_box(V);
      if (!preconditions(V, F, filename, controls))
        return;
    }
    spdlog::info("{} Pass Preconditions", filename);
    auto separate_type = PrismCage::SeparateType::kSurface;
    if (controls.at("allow-intersect"))
      separate_type = PrismCage::SeparateType::kNone;
    pc.reset(new PrismCage(V, F, 0.2, init_thick, separate_type));

    spdlog::enable_backtrace(42);
    if (!cage_is_away_from_ref(*pc)) {
      spdlog::error("initial intersect cage");
      exit(1);
    }
    spdlog::info("Initial Good.");
    spdlog::info("initialized, saving to {}.init", ser_file);
    pc->serialize(ser_file + ".init");
  }

  Vec3d bb = pc->ref.V.colwise().maxCoeff() - pc->ref.V.colwise().minCoeff();
  spdlog::info("Bounding Box {} {} {}, edge ratio {}", bb[0], bb[1], bb[2],
               edge_bb_ratio);
  double target_edge_length = bb.norm() * edge_bb_ratio;

  prism::local::RemeshOptions options(pc->mid.size(), target_edge_length);
  options.target_thickness = target_thickness;
  options.distortion_bound = angle_distortion;
  options.volume_centric = controls.at("volume-centric");
  options.parallel = controls.at("Parallel");

  if (controls.at("simplify-shell")) {
    remesh_schedule(*pc, options, ser_file);
  }
  spdlog::info("finished all passes, saving to {}", ser_file);

  if (controls.at("extract-stage")) {
    spdlog::info("saving to {}", ser_file);
    RowMatd mB, mT, V;
    RowMati F;
    vec2eigen(pc->mid, V);
    vec2eigen(pc->top, mT);
    vec2eigen(pc->base, mB);
    vec2eigen(pc->F, F);

    if (prism::cgal::polyhedron_self_intersect(mT, F)) {
      spdlog::warn("top {} self intersects", filename);
      prism::shell_extraction(*pc, false);
    }
    for (auto [v0, v1, v2] : pc->F) {
      if (v0 > v1 || v0 > v2) {
        spdlog::error("order error");
        exit(1);
      }
    }
    vec2eigen(pc->base, mB);
    vec2eigen(pc->F, F);
    if (prism::cgal::polyhedron_self_intersect(mB, F)) {
      spdlog::warn("base {} self intersects", filename);
      prism::shell_extraction(*pc, true);
    }
    for (auto [v0, v1, v2] : pc->F) {
      if (v0 > v1 || v0 > v2)
        spdlog::error("order error");
    }
  }
  pc->serialize(ser_file);
  spdlog::info("Memory Usage {}", getPeakRSS());
}

#include <igl/boundary_facets.h>
void boundary_perservation(std::string filename) {
  prism::geo::init_geogram();
  spdlog::info("Input H5 file {}", filename);
  PrismCage pc(filename);

  Vec3d bb = pc.ref.V.colwise().maxCoeff() - pc.ref.V.colwise().minCoeff();
  double target_edge_length = bb.norm();
  prism::local::RemeshOptions options(pc.mid.size(), target_edge_length);
  options.target_thickness = 1e-1;
  options.distortion_bound = 1e-3;

  std::map<std::pair<int, int>, std::vector<int>> meta_edges;
  std::vector<std::vector<int>> bnd;
  igl::boundary_loop(pc.ref.F, bnd); // this is assuming ref-bnd == mid-bnd
  for (auto &b : bnd) {
    for (int i = 0; i < b.size(); i++) {
      auto i1 = (i + 1) % b.size();
      meta_edges[std::pair(b[i], b[i + 1])] = {b[i], b[i + 1]};
    }
  }
  spdlog::set_level(spdlog::level::trace);

  // prism::local::featurecollapse_pass(pc, options, meta_edges);
  // prism::local::featurecollapse_pass(pc, options, meta_edges);
  spdlog::error("{}", meta_edges);
  pc.serialize("snail_save.h5");
}