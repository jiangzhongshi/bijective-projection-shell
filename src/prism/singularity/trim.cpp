#include <igl/boundary_loop.h>
#include <igl/copyleft/cgal/orient3D.h>
#include <igl/extract_manifold_patches.h>
#include <igl/per_face_normals.h>
#include <igl/read_triangle_mesh.h>
#include <igl/unique_edge_map.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/winding_number.h>
#include <igl/write_triangle_mesh.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>
#include <Eigen/Geometry>
#include <prism/cage_utils.hpp>
#include <prism/cgal/QP.hpp>
#include <prism/cgal/boolean.hpp>
#include <prism/common.hpp>
#include <prism/local_operations/local_mesh_edit.hpp>
#include <prism/local_operations/retain_triangle_adjacency.hpp>
#include <prism/local_operations/validity_checks.hpp>
#include <highfive/H5Easy.hpp>

RowMatd tetV =
    (RowMatd(4, 3) << 1, -1 / sqrt(3), -1 / sqrt(6), -1, -1 / sqrt(3),
     -1 / sqrt(6), 0, 2 / sqrt(3), -1 / sqrt(6), 0, 0, 3 / sqrt(6))
        .finished();
RowMati tetF = (RowMati(4, 3) << 0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2).finished();

namespace prism {
void remesh_intersect(const RowMatd& V1, const RowMati& F1, const RowMatd& V2,
                      const RowMati& F2, RowMatd& Vr, RowMati& Fr,
                      Eigen::VectorXi& CJ, Eigen::VectorXi& IM);

constexpr auto roll_shift_left = [](auto a, auto b, auto c,
                                    int s) -> std::tuple<int, int, int> {
  if (s == 0) return {a, b, c};
  if (s == 1)
    return {b, c, a};
  else
    return {c, a, b};
};

bool checks_orientation(const Vec3d& c, const Vec3d& p0, const Vec3d& p1,
                        const Vec3d& p2) {
  using igl::copyleft::cgal::orient3D;
  return orient3D(c.data(), p0.data(), p1.data(), p2.data()) > 0;
}

void extrude_intersection_curve(
    const Vec3d& center, const RowMatd& VN, const RowMatd& VC,
    const RowMati& FC, const prism::geogram::AABB& tree,
    const std::vector<std::vector<int>>& group_patches,  // rim, out
    const Eigen::VectorXi& intersection_curve, const RowMatd& cv_dirs,
    RowMatd& result, RowMati& rim_and_outer, Eigen::VectorXi& face_labels) {
  auto checks_for_triangle = [&center, &tree](const Vec3d& p0, const Vec3d& p1,
                                              const Vec3d& p2) -> bool {
    using igl::copyleft::cgal::orient3D;
    bool correct_orientation =
        orient3D(center.data(), p0.data(), p1.data(), p2.data()) > 0;

    bool intersection = false;  // tree.intersects_triangle({p0, p1, p2});
    return correct_orientation && !intersection;
  };
  auto rim = group_patches[0], out = group_patches[1];
  spdlog::info("Rim {} Out {}", rim.size(), out.size());
  std::set<int> fake_vertices;
  for (int i = 0; i < intersection_curve.size(); i++)
    fake_vertices.insert(intersection_curve[i]);
  for (auto f : out) {
    for (int j = 0; j < 3; j++) fake_vertices.insert(FC(f, j));
  }
  auto alpha = 1.;
  RowMatd outer_tip = VC;
  for (bool safe = false; !safe;) {
    safe = true;
    for (int i = 0; i < intersection_curve.size(); i++)
      outer_tip.row(intersection_curve[i]) =
          VC.row(intersection_curve[i]) + alpha * cv_dirs.row(i);
    for (auto i : out) {
      auto v0 = FC(i, 0), v1 = FC(i, 1), v2 = FC(i, 2);
      safe = checks_for_triangle(outer_tip.row(v0), outer_tip.row(v1),
                                 outer_tip.row(v2));
      if (!safe) {
        alpha *= 0.8;
        break;  // continue to next outer loop
      }
    }
  }
  spdlog::warn("alpha {}", alpha);
  // get a single alpha for all the cap,

  std::vector<double> face_steps(FC.rows(), 0.1);

  rim_and_outer.resize(rim.size() + out.size(),
                       3);  // face matrix for rim and out.
  {
    auto cnt = 0;
    for (auto& f : rim) rim_and_outer.row(cnt++) = FC.row(f);
    for (auto& f : out) rim_and_outer.row(cnt++) = FC.row(f);
  }
  // mark: number of fake vertices in each face: 1,2,3
  face_labels = Eigen::VectorXi::Ones(rim_and_outer.rows()) * 3;
  for (int f = 0; f < rim.size(); f++) {
    auto cnt = 0;
    std::vector<int> fake;
    for (int j = 0; j < 3; j++)
      if (fake_vertices.find(rim_and_outer(f, j)) != fake_vertices.end()) {
        fake.push_back(j);
      };
    face_labels[f] = fake.size();
    int shift = 0;  // shift_left, so that the fakes are in the front.
    if (fake.size() == 1) shift = fake[0];
    if (fake.size() == 2) {
      if (fake[1] - fake[0] == 1)
        shift = fake[0];
      else /*0,2*/
        shift = 2;
    };
    std::tie(rim_and_outer(f, 0), rim_and_outer(f, 1), rim_and_outer(f, 2)) =
        roll_shift_left(rim_and_outer(f, 0), rim_and_outer(f, 1),
                        rim_and_outer(f, 2), shift);
  }

  result = VC;
  for (auto v : fake_vertices) {
    result.row(v) = outer_tip.row(v);
  }
}

void corefine_with_tet(const RowMatd& mV, const RowMati& nbF, int v,
               /*output*/
               RowMatd& VC, RowMati& FC, Eigen::VectorXi& birth_face,
               std::vector<std::vector<int>>& group_patches,
               Eigen::VectorXi& intersection_curve,
               std::vector<int>& from_c_to_m) {
  group_patches.resize(4);
  Eigen::VectorXi one_ring_loop;
  igl::boundary_loop(nbF, one_ring_loop);
  double radius = 1;
  for (int i = 0; i < nbF.rows(); i++) {
    radius = std::min(radius, (mV.row(one_ring_loop[i]) - mV.row(v)).norm());
  }
  RowMatd local_tet = (tetV * radius).rowwise() + mV.row(v);
  for (auto tet_too_big = true; tet_too_big;) {
    radius /= 2;
    tet_too_big = false;
    local_tet = (tetV * radius).rowwise() + mV.row(v);
    // for (int j = 0; j < 4; j++) { TODO, need to modify interface for
    // one-ring
    //   if (tree.intersects_triangle({local_tet.row(tetF(j, 0)),
    //                                 local_tet.row(tetF(j, 1)),
    //                                 local_tet.row(tetF(j, 2))}),
    //                                 true) {
    //                                   spdlog::info("Inter");
    //     tet_too_big = true;
    //     break;
    //   }
    // }
  }

  spdlog::info("mV.rows() {}", mV.rows());
  spdlog::info("radius {}", radius);
  Eigen::VectorXi I;
  prism::remesh_intersect(mV, nbF, local_tet, tetF, VC, FC, birth_face, I);
  spdlog::info("VC.rows() {}", VC.rows());
  from_c_to_m.resize(VC.rows(), -1);
  for (int i = 0; i < VC.rows(); i++) {
    for (int j = 0; j < mV.rows(); j++)
      if (VC.row(i) == mV.row(j)) {  // use raw coordinate check.
        from_c_to_m[i] = j;
      }
  }

  {
    Eigen::VectorXi P;
    Eigen::MatrixXi E, uE;
    Eigen::VectorXi EMAP;
    std::vector<std::vector<int>> uE2E;
    igl::unique_edge_map(FC, E, uE, EMAP, uE2E);
    igl::extract_manifold_patches(FC, EMAP, uE2E, P);
    auto center_patch = -1, ring_rim = -1, tet_out = -1, tet_in = -1;
    for (int i = 0; i < P.rows(); i++) {
      assert(P[i] < 4);
      group_patches[P[i]].push_back(i);
      for (int j = 0; j < 3; j++) {
        if (VC.row(FC(i, j)) == mV.row(v))
          center_patch = P[i];
        else if (VC.row(FC(i, j)) == mV.row(one_ring_loop[0]))
          ring_rim = P[i];
      }
    }
    assert(center_patch != -1);
    assert(ring_rim != -1);
    // use the center faces, so the boundary loop is counterclockwise for
    // tet-out.

    RowMati center_faces(group_patches[center_patch].size(), 3);
    for (int i = 0; i < center_faces.rows(); i++)
      center_faces.row(i) = FC.row(group_patches[center_patch][i]);
    igl::boundary_loop(center_faces, intersection_curve);

    // erase center and put ring_rim first.
    group_patches.erase(group_patches.begin() + center_patch);
    if (center_patch < ring_rim) ring_rim--;
    if (ring_rim > 0) std::swap(group_patches[ring_rim], group_patches[0]);

    // then, differentiate between gp[1] and gp[2] to out and in.
    // edges in intersection_curve[0], intersection_curve[1];
    bool need_switch = false;
    auto ic0 = intersection_curve[0], ic1 = intersection_curve[1];
    for (auto f : group_patches[2]) {
      for (int j = 0; j < 3; j++) {
        if (ic0 == FC(f, j) && ic1 == FC(f, (j + 1) % 3)) need_switch = true;
      }
      if (need_switch == true) {  // curve is positive for 2-> out: need switch
        std::swap(group_patches[1], group_patches[2]);
        break;
      }
    }
  }
}


void trim(const RowMatd& mV, const RowMati& mF, const RowMatd& VN,
          const RowMatd& FN, const prism::geogram::AABB& tree, int v,
          const std::vector<int>& nb, RowMatd& resultV, RowMati& resultF,
          Eigen::VectorXi& C2M, Eigen::VectorXi& face_labels) {
  auto nb_num = nb.size();
  RowMatd nbFN(nb_num, 3);
  RowMati nbF(nb_num, 3);
  for (int i = 0; i < nb_num; i++) {
    nbFN.row(i) = FN.row(nb[i]);
    nbF.row(i) = mF.row(nb[i]);
  }

  spdlog::info("nbF,{}, v{}", nbF.rows(), v);
  RowMatd VC;  // C means combined, after boolean operation
  RowMati FC;
  Eigen::VectorXi birth_face;
  std::vector<std::vector<int>> group_patches(4);
  Eigen::VectorXi intersection_curve;
  std::vector<int> from_c_to_m;
  corefine_with_tet(mV, nbF, v,
       /*output*/
       VC, FC, birth_face, group_patches, intersection_curve, from_c_to_m);

  RowMatd VNC = RowMatd::Zero(VC.rows(), 3);
  for (int i = 0; i < VC.rows(); i++) VNC.row(i) = VN.row(from_c_to_m[i]);
  // valid directions for each point on the intersection curve (outwards)
  RowMatd curve_directions = RowMatd::Zero(intersection_curve.size(), 3);
  std::vector<std::vector<int>> cVF, cVFi;
  igl::vertex_triangle_adjacency(VC, FC, cVF, cVFi);
  for (int i = 0; i < intersection_curve.size(); i++) {
    auto cv = intersection_curve[i];
    Vec3d normal(0, 0, 0);
    std::set<int> faces_around;
    for (auto f : cVF[cv])
      if (birth_face[f] < nbF.rows()) faces_around.insert(birth_face[f]);
    assert(faces_around.size() <= 2);
    for (auto f : faces_around) normal += nbFN.row(f);
    curve_directions.row(i) = normal.stableNormalized();
  }

  // extrude intersection curve.
  extrude_intersection_curve(mV.row(v), VNC, VC, FC, tree, group_patches,
                             intersection_curve, curve_directions,
                             /*output:*/ resultV, resultF, face_labels);

  C2M = Eigen::Map<Eigen::VectorXi>(from_c_to_m.data(), from_c_to_m.size());
  Eigen::VectorXi loop_c;
  igl::boundary_loop(resultF, loop_c);
  double alpha = 1.0;  // back shrinking loop_c
  RowMatd tempV = resultV;
  for (auto clean = false; !clean;) {
    clean = true;
    for (int i = 0; i < loop_c.size(); i++) {
      assert(from_c_to_m[loop_c[i]] != -1 && "loop is wrong");
      resultV.row(loop_c[i]) =
          tempV.row(loop_c[i]) + alpha * VNC.row(loop_c[i]);
    }
    for (int i = 0; i < resultF.rows(); i++) {
      auto v0 = resultF(i, 0), v1 = resultF(i, 1), v2 = resultF(i, 2);
      clean = checks_orientation(mV.row(v), resultV.row(v0), resultV.row(v1),
                                 resultV.row(v2));
      if (!clean) {
        alpha *= 0.8;
        break;
      }
    }
  }
  spdlog::warn("loop alpha {}", alpha);
}

}  // namespace prism