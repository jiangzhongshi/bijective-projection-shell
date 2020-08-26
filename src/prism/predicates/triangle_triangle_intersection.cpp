/*
 *  Triangle-Triangle Overlap Test Routines
 *  July, 2002
 *  Updated December 2003
 *
 *  This file contains C implementation of algorithms for
 *  performing two and three-dimensional triangle-triangle intersection test
 *  The algorithms and underlying theory are described in
 *
 * "Fast and Robust Triangle-Triangle Overlap Test
 *  Using Orientation Predicates"  P. Guigue - O. Devillers
 *
 *  Journal of Graphics Tools, 8(1), 2003
 *
 *  Several geometric predicates are defined.  Their parameters are all
 *  points.  Each point is an array of two or three real precision
 *  floating point numbers. The geometric predicates implemented in
 *  this file are:
 *
 *
 *    int tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2,
 *                                     coplanar,source,target)
 *
 *       is a version that computes the segment of intersection when
 *       the triangles overlap (and are not coplanar)
 *
 *    each function returns 1 if the triangles (including their
 *    boundary) intersect, otherwise 0
 *
 *
 *  Other information are available from the Web page
 *  http:<i>//www.acm.org/jgt/papers/GuigueDevillers03/
 *
 */

// modified by Aaron to better detect coplanarity

#include "triangle_triangle_intersection.hpp"
typedef double real; // double

#include <geogram/numerics/predicates.h>
inline int sub_sub_cross_sub_dot(const real pa[3], const real pb[3],
                                 const real pc[3], const real pd[3]) {
  auto result = -GEO::PCK::orient_3d(pa, pb, pc, pd);
  if (result > 0)
    return 1;
  else if (result < 0)
    return -1;
  return 0;
}

inline int orient2d(const real a[3], const real b[3], const real c[3]) {
  auto result = GEO::PCK::orient_2d(a, b, c);
  // auto numer = ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]));
  // if (result >0) assert(numer>0);
  if (result > 0)
    return 1;
  else if (result < 0)
    return -1;
  return 0;
}
/* function prototype */

int tri_tri_intersection_test_3d(const real p1[3], const real q1[3],
                                 const real r1[3], const real p2[3],
                                 const real q2[3], const real r2[3],
                                 int *coplanar, real source[3], real target[3]);

int sub_sub_cross_sub_dot(const real a[3], const real b[3], const real c[3],
                          const real d[3]);

/* coplanar returns whether the triangles are coplanar
 *  source and target are the endpoints of the segment of
 *  intersection if it exists)
 */

/* some 3D macros */

#define CROSS(dest, v1, v2)                                                    \
  dest[0] = v1[1] * v2[2] - v1[2] * v2[1];                                     \
  dest[1] = v1[2] * v2[0] - v1[0] * v2[2];                                     \
  dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2)                                                      \
  dest[0] = v1[0] - v2[0];                                                     \
  dest[1] = v1[1] - v2[1];                                                     \
  dest[2] = v1[2] - v2[2];


/*
 *
 *  Two dimensional Triangle-Triangle Overlap Test
 *
 */

/* some 2D macros */

// #define orient2d(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))

#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2)                       \
  {                                                                            \
    if (orient2d(R2, P2, Q1) >= 0)                                          \
      if (orient2d(R2, Q2, Q1) <= 0)                                        \
        if (orient2d(P1, P2, Q1) > 0) {                                     \
          if (orient2d(P1, Q2, Q1) <= 0)                                    \
            return 1;                                                          \
          else                                                                 \
            return 0;                                                          \
        } else {                                                               \
          if (orient2d(P1, P2, R1) >= 0)                                    \
            if (orient2d(Q1, R1, P2) >= 0)                                  \
              return 1;                                                        \
            else                                                               \
              return 0;                                                        \
          else                                                                 \
            return 0;                                                          \
        }                                                                      \
      else if (orient2d(P1, Q2, Q1) <= 0)                                   \
        if (orient2d(R2, Q2, R1) <= 0)                                      \
          if (orient2d(Q1, R1, Q2) >= 0)                                    \
            return 1;                                                          \
          else                                                                 \
            return 0;                                                          \
        else                                                                   \
          return 0;                                                            \
      else                                                                     \
        return 0;                                                              \
    else if (orient2d(R2, P2, R1) >= 0)                                     \
      if (orient2d(Q1, R1, R2) >= 0)                                        \
        if (orient2d(P1, P2, R1) >= 0)                                      \
          return 1;                                                            \
        else                                                                   \
          return 0;                                                            \
      else if (orient2d(Q1, R1, Q2) >= 0) {                                 \
        if (orient2d(R2, R1, Q2) >= 0)                                      \
          return 1;                                                            \
        else                                                                   \
          return 0;                                                            \
      } else                                                                   \
        return 0;                                                              \
    else                                                                       \
      return 0;                                                                \
  };

#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2)                         \
  {                                                                            \
    if (orient2d(R2, P2, Q1) >= 0) {                                        \
      if (orient2d(P1, P2, Q1) >= 0) {                                      \
        if (orient2d(P1, Q1, R2) >= 0)                                      \
          return 1;                                                            \
        else                                                                   \
          return 0;                                                            \
      } else {                                                                 \
        if (orient2d(Q1, R1, P2) >= 0) {                                    \
          if (orient2d(R1, P1, P2) >= 0)                                    \
            return 1;                                                          \
          else                                                                 \
            return 0;                                                          \
        } else                                                                 \
          return 0;                                                            \
      }                                                                        \
    } else {                                                                   \
      if (orient2d(R2, P2, R1) >= 0) {                                      \
        if (orient2d(P1, P2, R1) >= 0) {                                    \
          if (orient2d(P1, R1, R2) >= 0)                                    \
            return 1;                                                          \
          else {                                                               \
            if (orient2d(Q1, R1, R2) >= 0)                                  \
              return 1;                                                        \
            else                                                               \
              return 0;                                                        \
          }                                                                    \
        } else                                                                 \
          return 0;                                                            \
      } else                                                                   \
        return 0;                                                              \
    }                                                                          \
  }

int ccw_tri_tri_intersection_2d(real p1[2], real q1[2], real r1[2], real p2[2],
                                real q2[2], real r2[2]) {
  if (orient2d(p2, q2, p1) >= 0) {
    if (orient2d(q2, r2, p1) >= 0) {
      if (orient2d(r2, p2, p1) >= 0)
        return 1;
      else
        INTERSECTION_TEST_EDGE(p1, q1, r1, p2, q2, r2)
    } else {
      if (orient2d(r2, p2, p1) >= 0)
        INTERSECTION_TEST_EDGE(p1, q1, r1, r2, p2, q2)
      else
        INTERSECTION_TEST_VERTEX(p1, q1, r1, p2, q2, r2)
    }
  } else {
    if (orient2d(q2, r2, p1) >= 0) {
      if (orient2d(r2, p2, p1) >= 0)
        INTERSECTION_TEST_EDGE(p1, q1, r1, q2, r2, p2)
      else
        INTERSECTION_TEST_VERTEX(p1, q1, r1, q2, r2, p2)
    } else
      INTERSECTION_TEST_VERTEX(p1, q1, r1, r2, p2, q2)
  }
};

int tri_tri_overlap_test_2d(real p1[2], real q1[2], real r1[2], real p2[2],
                            real q2[2], real r2[2]) {
  if (orient2d(p1, q1, r1) < 0)
    if (orient2d(p2, q2, r2) < 0)
      return ccw_tri_tri_intersection_2d(p1, r1, q1, p2, r2, q2);
    else
      return ccw_tri_tri_intersection_2d(p1, r1, q1, p2, q2, r2);
  else if (orient2d(p2, q2, r2) < 0)
    return ccw_tri_tri_intersection_2d(p1, q1, r1, p2, r2, q2);
  else
    return ccw_tri_tri_intersection_2d(p1, q1, r1, p2, q2, r2);
};

/*
 *
 *  Three-dimensional Triangle-Triangle Intersection
 *
 */

/*
   This macro is called when the triangles surely intersect
   It constructs the segment of intersection of the two triangles
   if they are not coplanar.
*/

constexpr auto contruct_intersection = [](auto p1, auto q1, auto r1, auto p2,
                                           auto q2, auto r2) {
  if (sub_sub_cross_sub_dot(q1, r2, p1, p2) > 0) {
    if (sub_sub_cross_sub_dot(r1, r2, p1, p2) <= 0) {
      if (sub_sub_cross_sub_dot(r1, q2, p1, p2) > 0) {
        return 1;
      } else {
        return 1;
      }
    } else {
      return 0;
    }
  } else {
    if (sub_sub_cross_sub_dot(q1, q2, p1, p2) < 0) {
      return 0;
    } else {
      if (sub_sub_cross_sub_dot(r1, q2, p1, p2) >= 0) {
        return 1;
      } else {
        return 1;
      }
    }
  }
};

/*
   The following version computes the segment of intersection of the
   two triangles if it exists.
   coplanar returns whether the triangles are coplanar
   source and target are the endpoints of the line segment of intersection
*/

// extern "C" real orient3d(const real *pa, const real *pb, const real *pc,
// const real *pd);

int coplanar_tri_tri3d(const real p1[3], const real q1[3], const real r1[3],
                       const real p2[3], const real q2[3], const real r2[3],
                       const real normal_1[3], const real normal_2[3]) {

  real P1[2], Q1[2], R1[2];
  real P2[2], Q2[2], R2[2];

  real n_x, n_y, n_z;

  n_x = ((normal_1[0] < 0) ? -normal_1[0] : normal_1[0]);
  n_y = ((normal_1[1] < 0) ? -normal_1[1] : normal_1[1]);
  n_z = ((normal_1[2] < 0) ? -normal_1[2] : normal_1[2]);

  /* Projection of the triangles in 3D onto 2D such that the area of
     the projection is maximized. */

  if ((n_x > n_z) && (n_x >= n_y)) {
    // Project onto plane YZ

    P1[0] = q1[2];
    P1[1] = q1[1];
    Q1[0] = p1[2];
    Q1[1] = p1[1];
    R1[0] = r1[2];
    R1[1] = r1[1];

    P2[0] = q2[2];
    P2[1] = q2[1];
    Q2[0] = p2[2];
    Q2[1] = p2[1];
    R2[0] = r2[2];
    R2[1] = r2[1];

  } else if ((n_y > n_z) && (n_y >= n_x)) {
    // Project onto plane XZ

    P1[0] = q1[0];
    P1[1] = q1[2];
    Q1[0] = p1[0];
    Q1[1] = p1[2];
    R1[0] = r1[0];
    R1[1] = r1[2];

    P2[0] = q2[0];
    P2[1] = q2[2];
    Q2[0] = p2[0];
    Q2[1] = p2[2];
    R2[0] = r2[0];
    R2[1] = r2[2];

  } else {
    // Project onto plane XY

    P1[0] = p1[0];
    P1[1] = p1[1];
    Q1[0] = q1[0];
    Q1[1] = q1[1];
    R1[0] = r1[0];
    R1[1] = r1[1];

    P2[0] = p2[0];
    P2[1] = p2[1];
    Q2[0] = q2[0];
    Q2[1] = q2[1];
    R2[0] = r2[0];
    R2[1] = r2[1];
  }

  return tri_tri_overlap_test_2d(P1, Q1, R1, P2, Q2, R2);
};

int tri_tri_intersection_test_3d(const real p1[3], const real q1[3],
                                 const real r1[3], const real p2[3],
                                 const real q2[3], const real r2[3],
                                 int *coplanar, real source[3], real target[3])

{
  int dp1, dq1, dr1, dp2, dq2, dr2;

  real v1[3], v2[3], v[3];
  real N1[3], N2[3], N[3];
  real alpha;

  SUB(v1, q1, p1)
  SUB(v2, r1, p1)
  CROSS(N1, v1, v2)

  SUB(v1, p2, r2)
  SUB(v2, q2, r2)
  CROSS(N2, v1, v2)

  *coplanar = 0;

  // Compute distance signs  of p1, q1 and r1
  // to the plane of triangle(p2,q2,r2)

  dp1 = sub_sub_cross_sub_dot(p2, q2, r2, p1);
  dq1 = sub_sub_cross_sub_dot(p2, q2, r2, q1);
  dr1 = sub_sub_cross_sub_dot(p2, q2, r2, r1);

  if (((dp1 * dq1) > 0) && ((dp1 * dr1) > 0))
    return 666;

  // Compute distance signs  of p2, q2 and r2
  // to the plane of triangle(p1,q1,r1)

  dp2 = sub_sub_cross_sub_dot(p1, q1, r1, p2);
  dq2 = sub_sub_cross_sub_dot(p1, q1, r1, q2);
  dr2 = sub_sub_cross_sub_dot(p1, q1, r1, r2);

  if (((dp2 * dq2) > 0) && ((dp2 * dr2) > 0))
    return 666;

  // Permutation in a canonical form of T1's vertices
  auto tri_tri_inter_3d =
      [&coplanar](const auto &p1, const auto &q1, const auto &r1,
                  const auto &p2, const auto &q2, const auto &r2,
                  const auto &dp2, const auto &dq2, const auto &dr2) {
        if (dp2 > 0) {
          if (dq2 > 0)
            return contruct_intersection(p1, r1, q1, r2, p2, q2);
          else if (dr2 > 0)
            return contruct_intersection(p1, r1, q1, q2, r2, p2);
          else
            return contruct_intersection(p1, q1, r1, p2, q2, r2);
        } else if (dp2 < 0) {
          if (dq2 < 0)
            return contruct_intersection(p1, q1, r1, r2, p2, q2);
          else if (dr2 < 0)
            return contruct_intersection(p1, q1, r1, q2, r2, p2);
          else
            return contruct_intersection(p1, r1, q1, p2, q2, r2);
        } else {
          if (dq2 < 0) {
            if (dr2 >= 0)
              return contruct_intersection(p1, r1, q1, q2, r2, p2);
            else
              return contruct_intersection(p1, q1, r1, p2, q2, r2);
          } else if (dq2 > 0) {
            if (dr2 > 0)
              return contruct_intersection(p1, r1, q1, p2, q2, r2);
            else
              return contruct_intersection(p1, q1, r1, q2, r2, p2);
          } else {
            if (dr2 > 0)
              return contruct_intersection(p1, q1, r1, r2, p2, q2);
            else if (dr2 < 0)
              return contruct_intersection(p1, r1, q1, r2, p2, q2);
            else {
              assert(false);
              *coplanar = 1;
              return -1;
            }
          }
        }
      };
  if (dp1 > 0) {
    if (dq1 > 0)
      return tri_tri_inter_3d(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2);
    else if (dr1 > 0)
      return tri_tri_inter_3d(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2);

    else
      return tri_tri_inter_3d(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2);
  } else if (dp1 < 0) {
    if (dq1 < 0)
      return tri_tri_inter_3d(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2);
    else if (dr1 < 0)
      return tri_tri_inter_3d(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2);
    else
      return tri_tri_inter_3d(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2);
  } else {
    if (dq1 < 0) {
      if (dr1 >= 0)
        return tri_tri_inter_3d(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2);
      else
        return tri_tri_inter_3d(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2);
    } else if (dq1 > 0) {
      if (dr1 > 0)
        return tri_tri_inter_3d(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2);
      else
        return tri_tri_inter_3d(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2);
    } else {
      if (dr1 > 0)
        return tri_tri_inter_3d(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2);
      else if (dr1 < 0)
        return tri_tri_inter_3d(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2);
      else {
        // triangles are co-planar

        *coplanar = 1;
        return coplanar_tri_tri3d(p1, q1, r1, p2, q2, r2, N1, N2);
      }
    }
  }
};

bool prism::predicates::triangle_triangle_overlap(
    const std::array<Vec3d, 3> &tri0, const std::array<Vec3d, 3> &tri1) {
  real src[3], trg[3];
  int cop[1];
  int flag = tri_tri_intersection_test_3d(
      tri0[0].data(), tri0[1].data(), tri0[2].data(), tri1[0].data(),
      tri1[1].data(), tri1[2].data(), cop, src, trg);
  if (flag == -1) {
    assert(false);
  };
  return flag == 1;
}

bool prism::predicates::segment_triangle_overlap(
    const std::array<Vec3d, 2> &seg, const std::array<Vec3d, 3> &tri) {
  // Warning: written in reference to segment_triangle_overlap in geogram.
  // To be tested.
  auto &[p, q] = seg;
  auto &[a, b, c] = tri;
  GEO::Sign abcp = GEO::PCK::orient_3d(a.data(), b.data(), c.data(), p.data());
  GEO::Sign abcq = GEO::PCK::orient_3d(a.data(), b.data(), c.data(), q.data());
  if (abcp == 0){ // project to 2d
  if(abcq == 0) {
    // coplanar_seg_tri_overlap()
  } else {
    // check if p inside abc
    // coplanar_point_in_tri(p, abc)
  }
  }
  if (abcq == 0) {
    // check if q inside abc
    // coplanar_point_in_tri(q, abc)
  }
  if (abcp == abcq)
    return false;
  GEO::Sign s1 = GEO::PCK::orient_3d(p.data(), q.data(), a.data(), b.data());
  GEO::Sign s2 = GEO::PCK::orient_3d(p.data(), q.data(), b.data(), c.data());
  if (s1 != s2)
    return false;
  GEO::Sign s3 = GEO::PCK::orient_3d(p.data(), q.data(), c.data(), a.data());
  return s2 == s3;
}