#ifndef PRISM_SINGULARITY_TRIM
#define PRISM_SINGULARITY_TRIM
#include "../geogram/AABB.hpp"
#include "../common.hpp"

namespace prism {
// Inputs:
//     mV: vertex positions
//     nbF: list of one-ring neighbor
//     v: id of center
// Outputs:
//     VC, FC: combined mesh (onering + tet)
//     birth_face: index from FC to [nbF, tet]
//     group_patches: ring-rim, tet-outer, tet-inner, ring-center
//     intersection_curve: id into VC, of ring-tet intersection curve
//     from_c_to_m: map VC to mV
void corefine_with_tet(const RowMatd& mV, const RowMati& nbF, int v,
                       /*output*/
                       RowMatd& VC, RowMati& FC, Eigen::VectorXi& birth_face,
                       std::vector<std::vector<int>>& group_patches,
                       Eigen::VectorXi& intersection_curve,
                       std::vector<int>& from_c_to_m);

// Inputs:
//    mV,mF, normals per vertex and face, id of singularity and its
//    neighbor
// Outputs:
//    resultV, resultF form a patch.
//    resultV contain all the boolean vertices.
//    resultF only rim+outercap, reordered based on fake vertices
//    C2M map from rV (VC) to mV,
//    face_labels: for each f in F, how many fake vertices are there. (1,2,3)
void trim(const RowMatd& mV, const RowMati& mF, const RowMatd& VN,
          const RowMatd& FN, const prism::geogram::AABB& tree, int v,
          const std::vector<int>& nb, RowMatd& resultV, RowMati& resultF,
          Eigen::VectorXi& C2M, Eigen::VectorXi& face_labels);
}  // namespace prism

#endif