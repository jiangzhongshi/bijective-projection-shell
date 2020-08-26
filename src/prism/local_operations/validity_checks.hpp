#ifndef PRISM_LOCAL_OPERATIONS_VALIDITY_CHECKS_HPP
#define PRISM_LOCAL_OPERATIONS_VALIDITY_CHECKS_HPP

#include "../geogram/AABB.hpp"
#include "../common.hpp"
#include "local_mesh_edit.hpp"
namespace prism{struct HashGrid;}
namespace prism::local_validity {

// this augments igl::volume numerical check with predicate check.
bool prism_positivity_with_numerical(const std::array<Vec3d, 6>& verts,
                              const std::array<bool, 3>& constrained = {
                                  false, false, false});

int attempt_relocate(
    std::vector<Vec3d>& base, std::vector<Vec3d>& top, std::vector<Vec3d>& mid,
    const std::vector<Vec3i>& F, const prism::geogram::AABB& tree,
    const RowMatd& refV, const RowMati& refF,
    const std::vector<std::set<int>>& map_track, double distortion_bound,
    // specified infos below
    const std::vector<int>& nb, int vid,
    const std::array<Vec3d, 3> /*b-m-t*/& relocations,
    std::vector<std::set<int>>& trackee);

int attempt_flip(const std::vector<Vec3d>& base, const std::vector<Vec3d>& top,
                 const std::vector<Vec3d>& mid, const std::vector<Vec3i>& F,
                 const prism::geogram::AABB& tree, const RowMatd& refV,
                 const RowMati& refF,
                 const std::vector<std::set<int>>& map_track,
                 double distortion_bound,
                 // specified infos below
                 int f0, int f1, int e0, int e1, int v0, int v1,
                 std::tuple<std::vector<int>,             /*shift*/
                            std::vector<std::set<int>> /*track*/
                            >&);

int attempt_split(
    std::vector<Vec3d>& base, std::vector<Vec3d>& top, std::vector<Vec3d>& mid,
    const std::vector<Vec3i>& F, const prism::geogram::AABB& tree,
    const RowMatd& refV, const RowMati& refF,
    const std::vector<std::set<int>>& map_track, double distortion_bound,
    bool improve_quality,
    // specified infos below
    int f0, int f1, int e0, int e1,
    const std::array<Vec3d, 3> /*b-m-t*/& newlocations, 
    std::tuple<std::vector<int> /*fid*/, std::vector<int>, /*shift*/
               std::vector<std::set<int>>               /*track*/
               >& checker);

int attempt_collapse(
    const std::vector<Vec3d>& base, const std::vector<Vec3d>& top,
    const std::vector<Vec3d>& mid, const std::vector<Vec3i>& F,
    const prism::geogram::AABB& tree, const std::array<std::shared_ptr<prism::HashGrid>,2> &grid,
    const RowMatd& refV, const RowMati& refF,
    const std::vector<std::set<int>>& map_track, double distortion_bound,
    double improve_quality_threshold,
    // specified infos below
    const std::vector<std::pair<int, int>>& neighbor0,
    const std::vector<std::pair<int, int>>& neighbor1, int f0, int f1, int u0,
    int u1,
    std::tuple<std::vector<int> /*newfid*/, std::vector<int> /*shifts*/,
               std::vector<std::set<int>> /*track*/>&);

// (1) if any new volume is negative
// requires to find vector<new Tri>, and their (modified) base-top V

bool volume_check(const std::vector<Vec3d>& base, const std::vector<Vec3d>& mid,
                  const std::vector<Vec3d>& top, const std::vector<Vec3i>& tris,
                  int num_cons = 0);
bool volume_check(const std::vector<Vec3d>& base,
                  const std::vector<Vec3d>& top, const std::vector<Vec3i>& tris,
                  int num_cons = 0);
// (2) if new prism intersect with ref-sheet
// same as (1)

bool intersect_check(const std::vector<Vec3d>& base,
                     const std::vector<Vec3d>& top,
                     const std::vector<Vec3i>& tris,
                     const prism::geogram::AABB& tree);

// (3) if distortion exceeds the bound, or flip
// from vector<new Tri>, need to call placement(),
// and update position, compute distortion for each Tri
}  // namespace prism::local_validity

#endif