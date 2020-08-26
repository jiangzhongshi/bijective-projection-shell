#include <igl/copyleft/cgal/wire_mesh.h>
#include <igl/readOBJ.h>
#include <igl/read_triangle_mesh.h>
#include <igl/volume.h>
#include <igl/write_triangle_mesh.h>
#include <spdlog/spdlog.h>

#include <doctest.h>

#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/local_operations/remesh_pass.hpp"
#include "prism/predicates/positive_prism_volume_12.hpp"

#define WIREPATH "../../nutshell_data/leg-intersect/"
TEST_CASE("wire mesh") {
  for (std::string filename : {"coarse_wire",  "ref_wire"}) {
  // for (std::string filename : {"qslim_wires", "shell_wires", "ref_wires"}) {
    RowMatd V;
    RowMati F;
    // this is because my special writing
    igl::read_triangle_mesh(WIREPATH + filename + ".obj", V, F);
    RowMati E = F.leftCols(2);
    RowMatd Vout;
    Eigen::VectorXi J;
    RowMati Fout;
    igl::copyleft::cgal::wire_mesh(V, E, 5e-3, 6, false, Vout, Fout, J);
    igl::write_triangle_mesh(WIREPATH + filename + "_mesh.obj", Vout, Fout);
  }
}