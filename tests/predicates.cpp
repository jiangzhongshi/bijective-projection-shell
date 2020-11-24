#include <doctest.h>
#include <igl/Timer.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <prism/cgal/triangle_triangle_intersection.hpp>
#include <prism/geogram/geogram_utils.hpp>

#include "prism/predicates/inside_octahedron.hpp"
#include "prism/predicates/inside_prism_tetra.hpp"
#include "prism/predicates/positive_prism_volume_12.hpp"
#include "prism/predicates/triangle_triangle_intersection.hpp"

TEST_CASE("tetra-sanity") {
  prism::geo::init_geogram();
  Vec3d t0 = STANDARD_PRISM[0];
  Vec3d t1 = STANDARD_PRISM[1];
  Vec3d t2 = STANDARD_PRISM[2];
  Vec3d t3 = STANDARD_PRISM[3];

  Vec3d ave = (t0 + t1 + t2 + t3) / 4;

  SUBCASE("Point in Tet") {
    // make sure the overall order is not screwed
    CHECK(prism::predicates::point_in_tetrahedron(ave, t0, t1, t2, t3));
  }

  SUBCASE("split type switch sanity") {
    // trivial type switch
    for (auto split : {true, false})
      CHECK(prism::predicates::point_in_prism(ave, split, STANDARD_PRISM));
  }

  SUBCASE("corners") {
    for (auto i : {0, 1, 2, 3, 4, 5}) {
      CAPTURE(i);
      CHECK(prism::predicates::point_in_prism(STANDARD_PRISM[i], false,
                                              STANDARD_PRISM));
    }
  }
}

TEST_CASE("volume-pos") {
  prism::geo::init_geogram();
  CHECK(prism::predicates::positive_prism_volume(STANDARD_PRISM));

  CHECK_FALSE(prism::predicates::positive_prism_volume(
      {STANDARD_PRISM[3], STANDARD_PRISM[4], STANDARD_PRISM[5],
       STANDARD_PRISM[0], STANDARD_PRISM[1], STANDARD_PRISM[2]}));
}

TEST_CASE("triangle intersect prism") {
  prism::geo::init_geogram();
  SUBCASE("Sanity, First Vertex is Center") {
    std::array<Vec3d, 3> tri;
    tri[0] = Vec3d(0, 0, 0);
    for (auto &a : STANDARD_PRISM) tri[0] += a / 6;
    CAPTURE(tri[0]);
    tri[1] = Vec3d(0, 1, 2);
    tri[2] = Vec3d(2, 5, 6);

    CHECK(prism::predicates::triangle_intersects_prism(tri, true,
                                                       STANDARD_PRISM));
  }

  SUBCASE("All vertices outside") {
    std::array<Vec3d, 3> tri;
    tri[0] = Vec3d(-1, -2, 0);
    tri[1] = Vec3d(2, 4, 0);
    tri[2] = Vec3d(-2, 4, 1);
    CAPTURE(tri[0]);

    CHECK(prism::predicates::triangle_intersects_prism(tri, true,
                                                       STANDARD_PRISM));
  }

  SUBCASE("Small All In") {
    std::array<Vec3d, 3> tri;
    tri[0] = Vec3d(2, 1, 0) / 100;
    tri[1] = Vec3d(2, 4, 0) / 100;
    tri[2] = Vec3d(2, 4, 1) / 100;
    CAPTURE(tri[0]);
    CAPTURE(tri[1]);
    CAPTURE(tri[2]);

    prism::geo::init_geogram();
    CHECK(prism::predicates::triangle_intersects_prism(tri, true,
                                                       STANDARD_PRISM));
  }
}

TEST_CASE("octahedron") {
  prism::geo::init_geogram();
  spdlog::set_level(spdlog::level::info);
  std::array<Vec3d, 3> base{STANDARD_PRISM[0], STANDARD_PRISM[1],
                            STANDARD_PRISM[2]};
  std::array<Vec3d, 3> top{STANDARD_PRISM[3], STANDARD_PRISM[4],
                           STANDARD_PRISM[5]};
  std::array<bool, 3> flag;
  prism::determine_convex_octahedron(base, top, flag, false);
  for (int i = 0; i < 3; i++) CHECK_FALSE(flag[i]);
  top[0][1] = -1;

  prism::determine_convex_octahedron(base, top, flag, false);
  CHECK(flag[0]);
  CHECK_FALSE(flag[1]);
  CHECK(flag[2]);
  {
    std::array<Vec3d, 3> tri;
    tri[0] = Vec3d(1, 2, 0) / 100;
    tri[1] = Vec3d(2, 4, 0) / 100;
    tri[2] = Vec3d(0.4, 0, 1.5);
    CAPTURE(tri[0]);

    CHECK(prism::triangle_intersect_octahedron(base, top, flag, tri, false));
  }
  {
    std::array<Vec3d, 3> tri;
    tri[0] = Vec3d(1, 2, 0);
    tri[1] = Vec3d(2, 4, 0);
    tri[2] = Vec3d(2, 4, 1);
    CAPTURE(tri[0]);

    CHECK_FALSE(
        prism::triangle_intersect_octahedron(base, top, flag, tri, false));

    //  CHECK_FALSE(prism::triangle_intersect_octahedron(base, top,  tri));
  }
}

TEST_CASE("octahedron_realdata") {
  prism::geo::init_geogram();
  RowMatd data(9, 3);
  data << -0.11822643141457159, -0.42141835414569256, 0.1349110937386112,
      0.015561042748273093, -0.7135878969126566, -0.13897739547590088,
      0.20534750976600258, -0.2746094405976865, 0.021047636295594176,
      -0.12054046929637627, -0.47763907937330585, 0.2846905158351154,
      0.1094766670461284, -0.7924777483949181, -0.0362337978853414,
      0.3164113699825023, -0.3173251291250236, 0.12800645504115316,
      -0.4539819319071452, -0.27586206896551724, 0.06797599176509482,
      -0.4539819319071452, -0.3103448275862069, 0.06797599176509482,
      -0.4344669454589205, -0.3085255729099737, 0.095271562358336;
  std::array<Vec3d, 3> base, top, tri;
  for (int i = 0; i < 3; i++) {
    base[i] = data.row(i);
    top[i] = data.row(i + 3);
    tri[i] = data.row(i + 6);
  }
  std::array<bool, 3> flag;
  prism::determine_convex_octahedron(base, top, flag, false);
  CHECK_FALSE(
      prism::triangle_intersect_octahedron(base, top, flag, tri, false));
}

TEST_CASE("tri-tri") {
  prism::geo::init_geogram();
  SUBCASE("equal tri") {
    std::array<Vec3d, 3> tri;
    tri[0] = Vec3d(-1, -2, 0);
    tri[1] = Vec3d(2, 4, 0);
    tri[2] = Vec3d(-2, 4, 1);
    auto tri2 = tri;
    CHECK(prism::predicates::triangle_triangle_overlap(tri, tri2));
    CHECK(prism::cgal::triangle_triangle_overlap(tri, tri2));
  }

  SUBCASE("corner") {
    std::array<Vec3d, 3> tri0{Vec3d(0.780465, -0.52344, -0.249586),
                              Vec3d(0, 0, 0),
                              Vec3d(-0.871657, 0.804416, 0.0250707)};
    std::array<Vec3d, 3> tri1{Vec3d(-0.959954, 0.70184, 0.335448),
                              Vec3d(0, 0, 0),
                              Vec3d(-0.873808, 0.0795207, -0.921439)};
    CHECK(prism::predicates::triangle_triangle_overlap(tri0, tri1));
    // CHECK(prism::cgal::triangle_triangle_overlap(tri0, tri1));
  }

  SUBCASE("coplanar") {
    std::array<Vec3d, 3> tri0{Vec3d(-1, -2, 0), Vec3d(2, 4, 0),
                              Vec3d(-2, 4, 0)};
    std::array<Vec3d, 3> tri1{Vec3d(10, 20, 0), Vec3d(3, 5, 0),
                              Vec3d(15, 40, 0)};
    CHECK_FALSE(prism::predicates::triangle_triangle_overlap(tri0, tri1));
    CHECK_FALSE(prism::cgal::triangle_triangle_overlap(tri0, tri1));
    std::array<Vec3d, 3> tri2{Vec3d(10, 20, 0), Vec3d(2, 4, 0),
                              Vec3d(15, 40, 0)};
    CHECK(prism::cgal::triangle_triangle_overlap(tri0, tri2));
    CHECK(prism::cgal::triangle_triangle_overlap(tri1, tri2));
    CHECK(prism::predicates::triangle_triangle_overlap(tri0, tri2));
    CHECK(prism::predicates::triangle_triangle_overlap(tri1, tri2));
  }
  SUBCASE("random") {
    Eigen::MatrixXd P(6, 3);
    for (int i = 0; i < 200; i++) {
      P.setRandom();
      if (i % 2 == 0) P.col(1).setZero();
      std::array<Vec3d, 3> tri0{P.row(0), P.row(1), P.row(2)};
      std::array<Vec3d, 3> tri1{P.row(3), P.row(4), P.row(5)};
      CHECK(prism::predicates::triangle_triangle_overlap(tri0, tri1) ==
            prism::cgal::triangle_triangle_overlap(tri0, tri1));
      CAPTURE(P);
      CHECK(prism::predicates::triangle_triangle_overlap(tri0, tri1) ==
            prism::cgal::triangle_triangle_overlap(tri0, tri1));
    }
  }
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersection_3.h>
bool segment_triangle_overlap(const std::array<Vec2d, 2> &seg,
                              const std::array<Vec2d, 3> &tri) {
  typedef ::CGAL::Exact_predicates_inexact_constructions_kernel K;
  std::array<K::Point_2, 5> cgal_points;
  for (int i = 0; i < 3; i++) cgal_points[i] = K::Point_2(tri[i][0], tri[i][1]);
  for (int i = 0; i < 2; i++)
    cgal_points[i + 3] = K::Point_2(seg[i][0], seg[i][1]);
  K::Triangle_2 t0(cgal_points[0], cgal_points[1], cgal_points[2]);
  K::Segment_2 s0(cgal_points[3], cgal_points[4]);
  return CGAL::do_intersect(t0, s0);
}

namespace coplanar {
extern std::array<Eigen::Matrix3i, 3> _global_record;
bool seg_tri_overlap(const std::tuple<Vec2d &, Vec2d &> &seg,
                     const std::tuple<Vec2d &, Vec2d &, Vec2d &> &tri);
}  // namespace coplanar
#include <geogram/numerics/predicates.h>
TEST_CASE("seg-tri-2d") {
  prism::geo::init_geogram();
  SUBCASE("sanity") {
    Eigen::MatrixXd P(6, 2);
    for (int i = 0; i < 10; i++) {
      P.setRandom();
      // if (i % 2 == 0) P.col(1).setZero();
      std::array<Vec2d, 2> tri0{P.row(0), P.row(1)};
      std::array<Vec2d, 3> tri1{P.row(0), P.row(1), P.row(5)};
      REQUIRE_EQ(coplanar::seg_tri_overlap({tri0[0], tri0[1]},
                                           {tri1[0], tri1[1], tri1[2]}),
                 true);
    }
  }
  SUBCASE("random") {
    Eigen::MatrixXd P(6, 2);
    for (int i = 0; i < 2000; i++) {
      P.setRandom();
      P = (P * 1e3).cast<int>().cast<double>();
      CAPTURE(i);
      CAPTURE(P);
      std::array<Vec2d, 2> tri0{P.row(0), P.row(1)};
      std::array<Vec2d, 3> tri1{P.row(3), P.row(5), P.row(4)};
      if (i % 2 == 0) tri1[0] = P.row(0) * 0.125 + P.row(1) * 0.875;
      if (i % 5 == 0) tri1[1] = P.row(0) * 0.5 + P.row(1) * 0.5;
      if (i % 3 == 0) tri1[2] = P.row(0) * 0.25 + P.row(1) * 0.75;
      if (GEO::PCK::orient_2d(tri1[0].data(), tri1[1].data(), tri1[2].data()) ==
          0)
        continue;
      REQUIRE_EQ(coplanar::seg_tri_overlap({tri0[0], tri0[1]},
                                           {tri1[0], tri1[1], tri1[2]}),
                 segment_triangle_overlap(tri0, tri1));
    }
    // spdlog::info("pqa<0 \n{}", coplanar::_global_record[0]);
    // spdlog::info("pqa=0 \n{}", coplanar::_global_record[1]);
    // spdlog::info("pqa>0 \n{}", coplanar::_global_record[2]);
  }
}

TEST_CASE("seg-tri") {
  prism::geo::init_geogram();
  SUBCASE("equal tri") {
    std::array<Vec3d, 3> tri;
    tri[0] = Vec3d(-1, -2, 0);
    tri[1] = Vec3d(2, 4, 0);
    tri[2] = Vec3d(-2, 4, 1);
    std::array<Vec3d, 2> tri2 = {tri[0], tri[1]};
    CHECK(prism::predicates::segment_triangle_overlap(tri2, tri));
    CHECK(prism::cgal::segment_triangle_overlap(tri2, tri));
  }

  SUBCASE("coplanar") {
    std::array<Vec3d, 2> tri0{Vec3d(-1, -2, 0), Vec3d(2, 4, 0)};
    std::array<Vec3d, 3> tri1{Vec3d(10, 20, 0), Vec3d(3, 5, 0),
                              Vec3d(15, 40, 0)};
    CHECK_FALSE(prism::predicates::segment_triangle_overlap(tri0, tri1));
    CHECK_FALSE(prism::cgal::segment_triangle_overlap(tri0, tri1));
    std::array<Vec3d, 3> tri2{Vec3d(10, 20, 0), Vec3d(2, 4, 0),
                              Vec3d(15, 40, 0)};
    CHECK(prism::cgal::segment_triangle_overlap(tri0, tri2));
    CHECK(prism::predicates::segment_triangle_overlap(tri0, tri2));
  }
  SUBCASE("random") {
    Eigen::MatrixXd P(6, 3);
    for (int i = 0; i < 200; i++) {
      P.setRandom();
      if (i % 2 == 0) P.col(1).setZero();
      std::array<Vec3d, 2> tri0{P.row(0), P.row(1)};
      std::array<Vec3d, 3> tri1{P.row(3), P.row(4), P.row(5)};
      if (i % 5 == 0) tri0[0] = (tri1[0] + tri1[1] + tri1[2] * 2) / 4;
      if (i % 7 == 0) std::swap(tri0[0], tri0[1]);
      CHECK(prism::predicates::segment_triangle_overlap(tri0, tri1) ==
            prism::cgal::segment_triangle_overlap(tri0, tri1));
      CAPTURE(P);
      CHECK(prism::predicates::segment_triangle_overlap(tri0, tri1) ==
            prism::cgal::segment_triangle_overlap(tri0, tri1));
    }
  }
}

TEST_CASE("tri-tri with seg") {
  prism::geo::init_geogram();
  constexpr auto tri_tri_test = [](const std::array<Vec3d, 3> &tri0,
                                   const std::array<Vec3d, 3> &tri1) {
    for (int i = 0; i < 3; i++) {
      if (prism::predicates::segment_triangle_overlap(
              {tri0[i], tri0[(i + 1) % 3]}, tri1))
        return true;
      if (prism::predicates::segment_triangle_overlap(
              {tri1[i], tri1[(i + 1) % 3]}, tri0))
        return true;
    }
    return false;
  };
  SUBCASE("equal tri") {
    std::array<Vec3d, 3> tri;
    tri[0] = Vec3d(-1, -2, 0);
    tri[1] = Vec3d(2, 4, 0);
    tri[2] = Vec3d(-2, 4, 1);
    auto tri2 = tri;
    CHECK(prism::predicates::triangle_triangle_overlap(tri, tri2));
    CHECK(tri_tri_test(tri, tri2));
  }

  SUBCASE("corner") {
    std::array<Vec3d, 3> tri0{Vec3d(0.780465, -0.52344, -0.249586),
                              Vec3d(0, 0, 0),
                              Vec3d(-0.871657, 0.804416, 0.0250707)};
    std::array<Vec3d, 3> tri1{Vec3d(-0.959954, 0.70184, 0.335448),
                              Vec3d(0, 0, 0),
                              Vec3d(-0.873808, 0.0795207, -0.921439)};
    CHECK(prism::predicates::triangle_triangle_overlap(tri0, tri1));
    // CHECK(tri_tri_test(tri0, tri1));
  }

  SUBCASE("coplanar") {
    std::array<Vec3d, 3> tri0{Vec3d(-1, -2, 0), Vec3d(2, 4, 0),
                              Vec3d(-2, 4, 0)};
    std::array<Vec3d, 3> tri1{Vec3d(10, 20, 0), Vec3d(3, 5, 0),
                              Vec3d(15, 40, 0)};
    CHECK_FALSE(prism::predicates::triangle_triangle_overlap(tri0, tri1));
    CHECK_FALSE(tri_tri_test(tri0, tri1));
    std::array<Vec3d, 3> tri2{Vec3d(10, 20, 0), Vec3d(2, 4, 0),
                              Vec3d(15, 40, 0)};
    CHECK(tri_tri_test(tri0, tri2));
    CHECK(tri_tri_test(tri1, tri2));
    CHECK(prism::predicates::triangle_triangle_overlap(tri0, tri2));
    CHECK(prism::predicates::triangle_triangle_overlap(tri1, tri2));
  }
  SUBCASE("random") {
    Eigen::MatrixXd P(6, 3);
    for (int i = 0; i < 200; i++) {
      P.setRandom();
      if (i % 2 == 0) P.col(1).setZero();
      std::array<Vec3d, 3> tri0{P.row(0), P.row(1), P.row(2)};
      std::array<Vec3d, 3> tri1{P.row(3), P.row(4), P.row(5)};
      CHECK(prism::predicates::triangle_triangle_overlap(tri0, tri1) ==
            tri_tri_test(tri0, tri1));
      CAPTURE(P);
      CHECK(prism::predicates::triangle_triangle_overlap(tri0, tri1) ==
            tri_tri_test(tri0, tri1));
    }
  }

  SUBCASE("timer") {
    Eigen::MatrixXd P(6, 3);
    double cnt = 0;
    spdlog::info("Timer 1");
    for (int i = 0; i < 2000; i++) {
      P.setRandom();
      if (i % 2 == 0) P.col(1).setZero();
      std::array<Vec3d, 3> tri0{P.row(0), P.row(1), P.row(2)};
      std::array<Vec3d, 3> tri1{P.row(3), P.row(4), P.row(5)};
      igl::Timer timer;
      timer.start();
      if (prism::predicates::triangle_triangle_overlap(tri0, tri1)) cnt++;
      cnt += timer.getElapsedTimeInMilliSec();
    }
    spdlog::info("Timer {}", cnt / 2000);
    cnt = 0;
    for (int i = 0; i < 2000; i++) {
      P.setRandom();
      if (i % 2 == 0) P.col(1).setZero();
      std::array<Vec3d, 3> tri0{P.row(0), P.row(1), P.row(2)};
      std::array<Vec3d, 3> tri1{P.row(3), P.row(4), P.row(5)};
      igl::Timer timer;
      timer.start();
      if (tri_tri_test(tri0, tri1)) cnt++;
      cnt += timer.getElapsedTimeInMilliSec();
    }
    spdlog::info("Timer {}", cnt / 2000);
  }
}

#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <highfive/H5Easy.hpp>
TEST_CASE("predicate debug") {
  prism::geo::init_geogram();
  std::string filename = "../build_clang/debug_predicate.h5";
  if (!std::filesystem::exists(filename)) return;
  H5Easy::File file(filename, H5Easy::File::ReadOnly);
  auto s = H5Easy::load<Vec3d>(file, "s");
  auto t = H5Easy::load<Vec3d>(file, "t");
  auto v0 = H5Easy::load<Vec3d>(file, "v0");
  auto v1 = H5Easy::load<Vec3d>(file, "v1");
  auto v2 = H5Easy::load<Vec3d>(file, "v2");
  CHECK(prism::cgal::segment_triangle_overlap({s, t}, {v0, v1, v2}));
  CHECK(prism::predicates::segment_triangle_overlap({s, t}, {v0, v1, v2}));
}
#include <prism/cgal/tetrahedron_tetrahedron_intersection.hpp>
#include <prism/predicates/tetrahedron_overlap.hpp>
TEST_CASE("tetrahedron predicate") {
  SUBCASE("timer") {
    Eigen::MatrixXd P(8, 3);
    double cnt = 0, tim = 0;
    spdlog::info("Timer 1");
    srand((unsigned int)(0));
    for (int i = 0; i < 2000; i++) {
      P.setRandom();
      if (i % 2 == 0) P.col(1).setZero();
      std::array<Vec3d, 4> tri0{P.row(0), P.row(1), P.row(2), P.row(6)};
      std::array<Vec3d, 4> tri1{P.row(3), P.row(4), P.row(5), P.row(7)};
      igl::Timer timer;
      timer.start();
      if (prism::predicates::tetrahedron_tetrahedron_overlap(tri0, tri1)) cnt++;
      // bool cgal_answer = prism::cgal::tetrahedron_tetrahedron_intersection(tri0, tri1);
      // REQUIRE(prism::predicates::tetrahedron_tetrahedron_overlap(tri0, tri1) ==cgal_answer);
      tim += timer.getElapsedTimeInMilliSec();
    }
    spdlog::info("Timer {}, C{}", tim, cnt / 2000);
    REQUIRE_EQ(cnt, 591);
    srand((unsigned int)(0));
    cnt = 0;
    tim = 0;
    for (int i = 0; i < 2000; i++) {
      P.setRandom();
      if (i % 2 == 0) P.col(1).setZero();
      std::array<Vec3d, 4> tri0{P.row(0), P.row(1), P.row(2), P.row(6)};
      std::array<Vec3d, 4> tri1{P.row(3), P.row(4), P.row(5), P.row(7)};
      igl::Timer timer;
      timer.start();
      if (prism::cgal::tetrahedron_tetrahedron_intersection(tri0, tri1)) cnt++;
      tim += timer.getElapsedTimeInMilliSec();
    }
    REQUIRE_EQ(cnt, 591);
    spdlog::info("Timer {}, C{}", tim, cnt / 2000);
  }
}