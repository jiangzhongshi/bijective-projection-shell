#include "positive_prism_volume_12.hpp"

#include <geogram/numerics/predicates.h>
#include <Eigen/Dense>

bool prism::predicates::positive_prism_volume(
    const std::array<Vec3d, 6> &verts) {
  using GEO::PCK::orient_3d;
  for (int i = 0; i < 12; i++) {
    if (orient_3d(verts[TWELVE_TETRAS[i][0]].data(),
                  verts[TWELVE_TETRAS[i][1]].data(),
                  verts[TWELVE_TETRAS[i][2]].data(),
                  verts[TWELVE_TETRAS[i][3]].data()) <= 0)
      return false;
  }
  return true;
}

bool prism::predicates::positive_prism_volume(
    const std::array<Vec3d, 6> &verts, const std::array<bool, 3> &constrained,
    bool numerical) {
  using GEO::PCK::orient_3d;
  for (int i = 0; i < 3; i++) {
    if (constrained[i])
      continue;
    for (int j = 0; j < 4; j++) {
      auto &tet = TWELVE_TETRAS[i * 4 + j];
      if (numerical) { // also check numerical validity
        RowMat3d local_verts;
        for (int k = 1; k < 4; k++)
          local_verts.row(k - 1) = verts[tet[k]] - verts[tet[0]];
        if (local_verts.determinant() <= 0)
          return false;
      }

      if (orient_3d(verts[tet[0]].data(), verts[tet[1]].data(),
                    verts[tet[2]].data(), verts[tet[3]].data()) <= 0) {
        return false;
      }
    }
  }
  return true;
}