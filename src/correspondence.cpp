#include <geogram/mesh/mesh_AABB.h>
#include <igl/exact_geodesic.h>
#include <igl/heat_geodesics.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/write_triangle_mesh.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>
#include <prism/cage_utils.hpp>
#include <prism/common.hpp>
#include <prism/geogram/geogram_utils.hpp>
#include <prism/phong/projection.hpp>
#include <highfive/H5Easy.hpp>
#include <igl/writePLY.h>
#include "prism/PrismCage.hpp"
#include "prism/geogram/AABB.hpp"
#include "prism/phong/query_correspondence.hpp"
#include <igl/Hit.h>

void transfer_pipeline(std::string shellfile, std::string section_file,
                       int upsample_level) {
  PrismCage pc(shellfile);

  H5Easy::File file(section_file, H5Easy::File::ReadWrite);
  auto qP = H5Easy::load<RowMatd>(file, "mV");
  auto pV = H5Easy::load<RowMatd>(file, "pV");
  auto pF = H5Easy::load<RowMati>(file, "pF");

  Eigen::VectorXi qfid;
  RowMatd quv;
  prism::correspond_bc(pc, pV, pF, qP, qfid, quv);

  H5Easy::dump(file, "fid", qfid, H5Easy::DumpMode::Overwrite);
  H5Easy::dump(file, "qUV", quv, H5Easy::DumpMode::Overwrite);
}