#pragma once

#include "../common.hpp"

struct PrismCage;
namespace prism::geogram {
struct AABB;
};
namespace igl {
struct Hit;
};
namespace prism {
void correspond_bc(const PrismCage &pc, const RowMatd &pxV, const RowMati &pxF,
                   const RowMatd &queryP, Eigen::VectorXi &queryF,
                   RowMatd &queryUV);
bool project_to_proxy_mesh(const std::array<Vec3d, 9> &stack,
                           const prism::geogram::AABB &pxtree, 
                           bool type, const Vec3d &spatial, igl::Hit &hit);
bool project_to_ref_mesh(const PrismCage &pc,
                         const std::vector<std::set<int>> &track_to_prism,
                         const std::vector<int> &tris, const Vec3d &point_value,
                         Vec3d &point_on_ref);
} // namespace prism