#ifndef PRISM_PREDICATES_INSIDE_OCTAHEDRON_HPP
#define PRISM_PREDICATES_INSIDE_OCTAHEDRON_HPP

#include <array>
#include "../common.hpp"
namespace prism {
bool inside_convex_octahedron(const std::array<Vec3d, 3>& base,
                              const std::array<Vec3d, 3>& top,
                              const Vec3d&);


void determine_convex_octahedron(const std::array<Vec3d, 3>& base,
                                 const std::array<Vec3d, 3>& top,
                                 std::array<bool, 3>& oct_type, bool degenerate=false);

// a (conservative test) of whether a triangle intersects a prism.
// (1) any tri-vert is contained inside
// otherwise, all vertices are outside
// (2) find out if triangle intersects with prism-boundary faces.
// (2') prismatic pillar intersects the triangle, or triangle edge intersects prismatic face.
bool triangle_intersect_octahedron(const std::array<Vec3d, 3>& base,
                              const std::array<Vec3d, 3>& top,
                              const std::array<bool, 3>& oct_type,
                              const std::array<Vec3d,3>& tri,
                              bool degenerate);

// similar to the previous, but ignore singularity.
// Note: this is actually triangle vs. pyramid
bool singularless_triangle_intersect_octahedron(const std::array<Vec3d, 3>& base,
                              const std::array<Vec3d, 3>& top,
                              const std::array<bool, 3>& oct_type,
                              const std::array<Vec3d,3>& tri);

// intersection predicate, ignoring the first point
bool pointless_triangle_intersect_octahedron(const std::array<Vec3d, 3>& base,
                              const std::array<Vec3d, 3>& top,
                              const std::array<bool, 3>& oct_type,
                              const std::array<Vec3d,3>& tri);
}
#endif