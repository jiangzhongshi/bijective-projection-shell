
#include "AABB.hpp"

#include <geogram/basic/geometry_nd.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/numerics/predicates.h>
#include <igl/Hit.h>
#include <igl/Timer.h>
#include <igl/boundary_facets.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/squared_distance_3.h>

#include <iostream>
#include <limits>
#include <list>
#include <prism/geogram/geogram_utils.hpp>

#include "prism/predicates/triangle_triangle_intersection.hpp"

prism::geogram::AABB::AABB(const RowMatd &V, const RowMati &F, bool _enabled)
    : enabled(_enabled) {
  if (!enabled)
    return;

  geo_polyhedron_ptr_ = std::make_unique<GEO::Mesh>();
  prism::geo::to_geogram_mesh(V, F, *geo_polyhedron_ptr_);
  geo_tree_ptr_ =
      std::make_unique<GEO::MeshFacetsAABB>(*geo_polyhedron_ptr_, true);

  geo_vertex_ind.resize(V.rows());
  GEO::Attribute<int> original_indices(
      geo_polyhedron_ptr_->vertices.attributes(), "vertex_id");
  for (int i = 0; i < original_indices.size(); i++)
    geo_vertex_ind[i] = original_indices[i];

  geo_face_ind.resize(F.rows());
  GEO::Attribute<int> face_indices(geo_polyhedron_ptr_->facets.attributes(),
                                   "facet_id");
  for (int i = 0; i < face_indices.size(); i++)
    geo_face_ind[i] = face_indices[i];
}

bool prism::geogram::AABB::intersects_triangle(const std::array<Vec3d, 3> &P,
                                               bool use_freeze) const {
  using K = ::CGAL::Exact_predicates_inexact_constructions_kernel;
  if (!enabled)
    return false;
  using namespace GEO;
  Box in_box;
  for (int i = 0; i < 3; i++) {
    in_box.xyz_max[i] = std::max({P[0][i], P[1][i], P[2][i]});
    in_box.xyz_min[i] = std::min({P[0][i], P[1][i], P[2][i]});
  }

  K::Point_3 kp0(P[0][0], P[0][1], P[0][2]);
  K::Point_3 kp1(P[1][0], P[1][1], P[1][2]);
  K::Point_3 kp2(P[2][0], P[2][1], P[2][2]);
  K::Triangle_3 kquery(kp0, kp1, kp2);
  auto &ind = geo_vertex_ind;
  bool already_found = false;

  auto action = [&P, &pp = geo_polyhedron_ptr_, &already_found, &ind, &kquery,
                 &kp0, &kp1, &kp2, num_freeze = num_freeze,
                 use_freeze](GEO::index_t f) {
    if (already_found)
      return;

    index_t c = pp->facets.corners_begin(f);
    auto i0 = pp->facet_corners.vertex(c), i1 = pp->facet_corners.vertex(c + 1),
         i2 = pp->facet_corners.vertex(c + 2);
    assert(ind[i0] < ind[i1] && ind[i0] < ind[i2]);
    const vec3 &v0 = Geom::mesh_vertex(*pp, i0);
    const vec3 &v1 = Geom::mesh_vertex(*pp, i1);
    const vec3 &v2 = Geom::mesh_vertex(*pp, i2);
    if (use_freeze && ind[i0] < num_freeze) { // triangle vs segment
      K::Point_3 kv0(v0.x, v0.y, v0.z);
      K::Point_3 kv1(v1.x, v1.y, v1.z);
      K::Point_3 kv2(v2.x, v2.y, v2.z);
      if (kv1 == kp0)
        std::swap(kv1, kv0);
      if (kv2 == kp0)
        std::swap(kv2, kv0);
      assert(kv0 == kp0); // floating point should be exact here.
      K::Segment_3 ks(kv1, kv2);
      K::Triangle_3 kt(kv0, kv1, kv2);
      K::Segment_3 ks1(kp1, kp2);

      spdlog::trace("testing {} {} {}", ind[i0], ind[i1], ind[i2]);

      if (CGAL::do_intersect(ks, kquery) || CGAL::do_intersect(kt, ks1)) {
        already_found = true;
      }
      return;
    } // ordinary check triangle vs triangle
    Vec3d kv0(v0.x, v0.y, v0.z);
    Vec3d kv1(v1.x, v1.y, v1.z);
    Vec3d kv2(v2.x, v2.y, v2.z);
    std::array<Vec3d, 3> kt{kv0, kv1, kv2};
    if (prism::predicates::triangle_triangle_overlap(kt, P)) {
      spdlog::trace("kt {} {} {}",
                    kv0.format(Eigen::IOFormat(Eigen::FullPrecision)),
                    kv1.format(Eigen::IOFormat(Eigen::FullPrecision)),
                    kv2.format(Eigen::IOFormat(Eigen::FullPrecision)));
      already_found = true;
    }
  };
  geo_tree_ptr_->compute_bbox_facet_bbox_intersections(in_box, action);
  return already_found;
}

bool prism::geogram::AABB::segment_query(const Vec3d &start, const Vec3d &end,
                                         int &face_id,
                                         Vec3d &finalpoint) const {
  using K = ::CGAL::Exact_predicates_inexact_constructions_kernel;
  if (!enabled)
    return false;
  // returns the geogram index after reordering, and one need to map it back.
  K::Point_3 s(start[0], start[1], start[2]);
  K::Point_3 e(end[0], end[1], end[2]);
  K::Segment_3 cgal_seg(s, e);
  using namespace GEO;
  Box in_box;
  for (int i = 0; i < 3; i++) {
    in_box.xyz_max[i] = std::max({start[i], end[i]});
    in_box.xyz_min[i] = std::min({start[i], end[i]});
  }
  auto &pp = geo_polyhedron_ptr_;
  auto action = [&finalpoint, &start, &end, &cgal_seg, &pp,
                 &face_id](GEO::index_t f) {
    if (face_id >= 0)
      return;
    index_t c = pp->facets.corners_begin(f);
    const vec3 &p0 = Geom::mesh_vertex(*pp, pp->facet_corners.vertex(c));
    const vec3 &p1 = Geom::mesh_vertex(*pp, pp->facet_corners.vertex(c + 1));
    const vec3 &p2 = Geom::mesh_vertex(*pp, pp->facet_corners.vertex(c + 2));
    Vec3d kp0(p0.x, p0.y, p0.z);
    Vec3d kp1(p1.x, p1.y, p1.z);
    Vec3d kp2(p2.x, p2.y, p2.z);
    std::array<Vec3d, 3> kt{kp0, kp1, kp2};
    K::Point_3 cp0(p0.x, p0.y, p0.z);
    K::Point_3 cp1(p1.x, p1.y, p1.z);
    K::Point_3 cp2(p2.x, p2.y, p2.z);
    K::Triangle_3 c_kt(cp0, cp1, cp2);
    auto inter = CGAL::intersection(c_kt, cgal_seg);
    if (inter) {
      spdlog::trace("kt {} {} {}",
                    kp0.format(Eigen::IOFormat(Eigen::FullPrecision)),
                    kp1.format(Eigen::IOFormat(Eigen::FullPrecision)),
                    kp2.format(Eigen::IOFormat(Eigen::FullPrecision)));
      const K::Point_3 *point;
      if (inter.value().which() == 0) {
        point = &boost::get<K::Point_3>((inter).value());
      } else {
        point = &boost::get<K::Segment_3>((inter).value()).start();
      }
      finalpoint =
          Vec3d(CGAL::to_double(point->x()), CGAL::to_double(point->y()),
                CGAL::to_double(point->z()));
      face_id = f;
    }
  };
  geo_tree_ptr_->compute_bbox_facet_bbox_intersections(in_box, action);
  return face_id >= 0;
}

std::optional<Vec3d>
prism::geogram::AABB::segment_query(const Vec3d &start,
                                    const Vec3d &end) const {
  int face_id = -1;
  Vec3d finalpoint;
  segment_query(start, end, face_id, finalpoint);
  if (face_id != -1)
    return finalpoint;
  return {};
}

bool prism::geogram::AABB::segment_hit(const Vec3d &start, const Vec3d &end,
                                       igl::Hit &hit) const {
  int fid = -1;
  Vec3d finalpoint;
  if (!segment_query(start, end, fid, finalpoint))
    return false;
  auto &pp = geo_polyhedron_ptr_;
  using namespace GEO;
  index_t c = pp->facets.corners_begin(fid);
  const vec3 &p0 = Geom::mesh_vertex(*pp, pp->facet_corners.vertex(c));
  const vec3 &p1 = Geom::mesh_vertex(*pp, pp->facet_corners.vertex(c + 1));
  const vec3 &p2 = Geom::mesh_vertex(*pp, pp->facet_corners.vertex(c + 2));
  vec3 gp(finalpoint[0], finalpoint[1], finalpoint[2]);
  auto a0 = Geom::triangle_area(gp, p1, p2);
  auto a1 = Geom::triangle_area(gp, p2, p0);
  auto a2 = Geom::triangle_area(gp, p0, p1);
  assert(abs(Geom::triangle_area(p0, p1, p2) - (a0 + a1 + a2)) < 1e-10);
  hit.u = a1 / (a0 + a1 + a2);
  hit.v = a2 / (a0 + a1 + a2);
  hit.id = geo_face_ind[fid];
  hit.gid = 0;
  return true;
}

double prism::geogram::AABB::ray_length(const Vec3d &start, const Vec3d &dir,
                                        double max_step, int ignore_v) const {
  using K = CGAL::Exact_predicates_exact_constructions_kernel;
  if (!enabled)
    return max_step;
  // returns the geogram index after reordering, and one need to map it back.
  auto end = start + max_step * dir;
  K::Point_3 s(start[0], start[1], start[2]);
  K::Point_3 e(end[0], end[1], end[2]);
  K::Segment_3 cgal_seg(s, e);
  GEO::Box in_box;
  for (int i = 0; i < 3; i++) {
    in_box.xyz_max[i] = std::max({start[i], end[i]});
    in_box.xyz_min[i] = std::min({start[i], end[i]});
  }
  auto &pp = geo_polyhedron_ptr_;
  double min_length_sq = max_step * max_step;
  spdlog::trace("input {}", s);
  auto action = [&ignore_v, &s, &start, &end, &cgal_seg, &pp,
                 &geo_vertex_ind = geo_vertex_ind,
                 &min_length_sq](GEO::index_t f) -> void {
    using namespace GEO;
    index_t c = pp->facets.corners_begin(f);
    auto i0 = pp->facet_corners.vertex(c), i1 = pp->facet_corners.vertex(c + 1),
         i2 = pp->facet_corners.vertex(c + 2);
    if (geo_vertex_ind[i0] == ignore_v || geo_vertex_ind[i1] == ignore_v ||
        geo_vertex_ind[i2] == ignore_v) { // ignore itself.
      spdlog::trace("ignored.");
      return;
    }
    const vec3 &p0 = Geom::mesh_vertex(*pp, i0);
    const vec3 &p1 = Geom::mesh_vertex(*pp, i1);
    const vec3 &p2 = Geom::mesh_vertex(*pp, i2);
    Vec3d kp0(p0.x, p0.y, p0.z);
    Vec3d kp1(p1.x, p1.y, p1.z);
    Vec3d kp2(p2.x, p2.y, p2.z);
    std::array<Vec3d, 3> kt{kp0, kp1, kp2};
    K::Point_3 cp0(p0.x, p0.y, p0.z);
    K::Point_3 cp1(p1.x, p1.y, p1.z);
    K::Point_3 cp2(p2.x, p2.y, p2.z);
    K::Triangle_3 c_kt(cp0, cp1, cp2);
    auto inter = CGAL::intersection(c_kt, cgal_seg);
    if (inter) {
      spdlog::trace("kt {} {} {}",
                    kp0.format(Eigen::IOFormat(Eigen::FullPrecision)),
                    kp1.format(Eigen::IOFormat(Eigen::FullPrecision)),
                    kp2.format(Eigen::IOFormat(Eigen::FullPrecision)));
      if (inter.value().which() == 1)
        return; // parallel, another triangle will see this.
      const K::Point_3 *point = &boost::get<K::Point_3>((inter).value());
      min_length_sq = std::min(
          CGAL::to_double(CGAL::squared_distance(*point, s)), min_length_sq);
      spdlog::trace("point {}", *point);
      spdlog::trace("len {}", min_length_sq);
    }
  };
  geo_tree_ptr_->compute_bbox_facet_bbox_intersections(in_box, action);
  return std::sqrt(min_length_sq);
}

bool prism::geogram::AABB::numerical_self_intersection(double tol) const {
  // this function is oblivious to re-ordering (unless debugging), so no need to
  // include permutation in the regular code
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using namespace GEO;
  auto &pp = geo_polyhedron_ptr_;
  for (auto v = 0; v < geo_vertex_ind.size(); v++) {
    auto p = GEO::Geom::mesh_vertex(*pp, v);
    GEO::Box in_box;
    for (int i = 0; i < 3; i++) {
      in_box.xyz_max[i] = p[i] + tol;
      in_box.xyz_min[i] = p[i] - tol;
    }
    bool intersect_flag = false;
    K::Point_3 kp(p.x, p.y, p.z);
    auto action = [&pp, &p, &v, &intersect_flag, &kp,
                   &tol](GEO::index_t f) -> void {
      if (intersect_flag)
        return;
      GEO::index_t c = pp->facets.corners_begin(f);
      auto i0 = pp->facet_corners.vertex(c),
           i1 = pp->facet_corners.vertex(c + 1),
           i2 = pp->facet_corners.vertex(c + 2);
      // spdlog::trace("io {} {} {} v {}",i0,i1,i2,v);
      if (i0 == v || i1 == v || i2 == v) // is neighbor
        return;
      const vec3 &v0 = Geom::mesh_vertex(*pp, i0);
      const vec3 &v1 = Geom::mesh_vertex(*pp, i1);
      const vec3 &v2 = Geom::mesh_vertex(*pp, i2);
      K::Point_3 kv0(v0.x, v0.y, v0.z);
      K::Point_3 kv1(v1.x, v1.y, v1.z);
      K::Point_3 kv2(v2.x, v2.y, v2.z);
      K::Triangle_3 kt(kv0, kv1, kv2);
      if (CGAL::squared_distance(kp, kt) <= tol * tol) {
        intersect_flag = true;
        return;
      };
      return;
    };
    geo_tree_ptr_->compute_bbox_facet_bbox_intersections(in_box, action);
    if (intersect_flag)
      return true;
  }
  return false;
}
