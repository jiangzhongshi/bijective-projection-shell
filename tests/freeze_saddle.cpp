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
#include <prism/PrismCage.hpp>
#include <prism/cage_utils.hpp>
#include <prism/cgal/QP.hpp>
#include <prism/cgal/boolean.hpp>
#include <prism/extraction.hpp>
#include <prism/local_operations/local_mesh_edit.hpp>
#include <prism/local_operations/retain_triangle_adjacency.hpp>
#include <prism/singularity/trim.hpp>
#include <doctest.h>
#include <highfive/H5Easy.hpp>
#include "test_common.hpp"

TEST_CASE("trim shell") {
  RowMatd mV;
  RowMati mF;
  igl::read_triangle_mesh("/home/zhongshi/1146172_sf.obj", mV, mF);
  put_in_unit_box(mV);
  PrismCage pc(mV, mF, 0.2);
  for (int i = 0; i < pc.mid.size(); i++) {
    pc.base[i] = (pc.mid[i] + pc.base[i]) / 2;
    pc.top[i] = (pc.mid[i] + pc.top[i]) / 2;
  }
  spdlog::set_level(spdlog::level::debug);
  prism::shell_extraction(pc, false);
  prism::shell_extraction(pc, true);
  pc.serialize("extract.h5");
}

TEST_CASE("continue trim top shell") {
  PrismCage pc("extract.h5");
  pc.ref.aabb->num_freeze = 29;
  for (int i = 0; i < pc.mid.size(); i++) {
    while ((pc.base[i]-pc.mid[i]).norm() > 1e-2)
      pc.base[i] = 0.2*pc.mid[i] + 0.8*(pc.base[i]);
    while ((pc.top[i]-pc.mid[i]).norm() > 1e-2)
      pc.top[i] = 0.2*pc.mid[i] + 0.8*(pc.top[i]);
  }
  RowMatd mV;
  RowMati mF;
  vec2eigen(pc.top, mV);
  vec2eigen(pc.F, mF);
  std::vector<std::vector<int>> VF, VFi;
  igl::vertex_triangle_adjacency(mV, mF, VF, VFi);
  std::vector<Vec3d> finalV;
  std::vector<Vec3i> finalF;
  using Vec4i = std::array<int, 4>;
  std::vector<Vec4i> finalT;
  eigen2vec(mV, finalV);
  eigen2vec(mF, finalF);
  std::vector<int> birthV;
  for (int i = 0; i < mV.rows(); i++) birthV.push_back(i);

  for (int i = 0; i < pc.ref.aabb->num_freeze; i++) {
    RowMati nbF(VF[i].size(), 3);
    for (int j = 0; j < VF[i].size(); j++) {
      nbF.row(j) = mF.row(VF[i][j]);
    }
    RowMatd VC;
    RowMati FC;
    Eigen::VectorXi birthface, intersection_curve;
    std::vector<std::vector<int>> group_patches;
    std::vector<int> from_c_to_m;
    prism::corefine_with_tet(mV, nbF, i, VC, FC, birthface, group_patches,
                             intersection_curve, from_c_to_m);

    for (int c = 0; c < from_c_to_m.size(); c++) {
      auto& m = from_c_to_m[c];
      if (m == -1) {
        m = finalV.size();
        finalV.push_back(VC.row(c));
        birthV.push_back(i);
      }
    }
    std::vector<int> rim_outer;
    std::merge(group_patches[0].begin(), group_patches[0].end(),
    group_patches[1].begin(), group_patches[1].end(),std::inserter(rim_outer, rim_outer.end()));
    for (auto f : rim_outer) {
      int n0 = from_c_to_m[FC(f, 0)], n1 = from_c_to_m[FC(f, 1)],
          n2 = from_c_to_m[FC(f, 2)];
      finalF.emplace_back(Vec3i{n0, n1, n2});
      // finalT.emplace_back(Vec4i{i, n0, n1, n2});
    }
  }

  {
    RowMati fT;
    RowMatd fV;
    RowMati fF;
    vec2eigen(finalV, fV);
    vec2eigen(finalF, fF);
    vec2eigen(finalT, fT);
    H5Easy::File file("trim_top_shell.h5", H5Easy::File::Overwrite);
    H5Easy::dump(file, "mV", mV);
    H5Easy::dump(file, "mF", mF);
    H5Easy::dump(file, "fV", fV);
    H5Easy::dump(file, "fF", fF);
    H5Easy::dump(file, "fT", fT);
    H5Easy::dump(file, "birth", birthV);
  }
}

TEST_CASE("trim tet at singularity") {
  spdlog::set_level(spdlog::level::debug);
  RowMatd mV;
  RowMati mF;
  igl::read_triangle_mesh("/home/zhongshi/1146172_sf.obj", mV, mF);
  put_in_unit_box(mV);
  // igl::read_triangle_mesh("../tests/data/saddle/original.obj", mV, mF);
  std::vector<Vec3d> V;
  std::vector<Vec3i> F;
  eigen2vec(mV, V);
  eigen2vec(mF, F);
  auto [FF, FFi] = prism::local::triangle_triangle_adjacency(F);
  std::set<int> omni_singu;
  std::set<int> bnd;
  RowMatd VN, FN;
  bool good_normal =
      prism::cage_utils::most_visible_normals(mV, mF, VN, omni_singu);
  prism::cage_utils::reorder_singularity_to_front(mV, mF, VN, omni_singu,bnd);
  for (int i = 0; i < mF.rows(); i++) {  // reorder F
    auto [s, mt, shift] = tetra_split_AorB({mF(i, 0), mF(i, 1), mF(i, 2)});
    for (int j = 0; j < 3; j++) mF(i, j) = mt[j];
  }

  prism::geogram::AABB tree(mV, mF);
  tree.num_freeze = omni_singu.size();
  std::vector<double> face_steps = prism::cage_utils::extrude_along_normals(
      mV, mF, tree, VN, true, omni_singu.size(), 1e-4);
  std::vector<double> vertex_steps(mV.rows(), 1);
  for (int i = 0; i < mF.rows(); i++)
    for (int j = 0; j < 3; j++)
      vertex_steps[mF(i, j)] = std::min(vertex_steps[mF(i, j)], face_steps[i]);
  prism::cage_utils::iterative_retract_normal(mV, mF, tree, VN, true,
                                              omni_singu.size(), vertex_steps);
  igl::per_face_normals_stable(mV, mF, FN);

  spdlog::set_level(spdlog::level::trace);
  CHECK_FALSE(good_normal);

  std::vector<std::vector<int>> VF, VFi;
  igl::vertex_triangle_adjacency(mV, mF, VF, VFi);
  std::vector<Vec3d> finalV;
  std::vector<Vec3i> finalF;
  using Vec4i = std::array<int, 4>;
  std::vector<Vec4i> finalT;
  eigen2vec(mV, finalV);
  eigen2vec(mF, finalF);

  using igl::copyleft::cgal::orient3D;
  auto checks_orientation = [](const Vec3d& c, const Vec3d& p0, const Vec3d& p1,
                               const Vec3d& p2) -> bool {
    return orient3D(c.data(), p0.data(), p1.data(), p2.data()) > 0;
  };

  std::set<int, std::greater<int>> faces_to_erase;
  std::vector<int> birthV;
  for (int i = 0; i < mV.rows(); i++) birthV.push_back(i);

  for (int v = 0; v < omni_singu.size(); v++) {
    spdlog::warn("v: {}", v);
    RowMatd rV;
    RowMati rF;  // rim and out patch.
    Eigen::VectorXi rV2mV, rF_label;
    prism::trim(mV, mF, VN, FN, tree, v, VF[v], rV, rF, rV2mV, rF_label);
    for (int i = 0; i < rV2mV.rows(); i++) {
      auto mi = rV2mV[i];
      if (mi != -1) {  // old vert, need to compare
        if (rV.row(i) != mV.row(mi)) {
          double step = (rV.row(i) - mV.row(mi)).norm() / VN.row(mi).norm();
          vertex_steps[mi] = std::min(vertex_steps[mi], step);
        };
        finalV[mi] = mV.row(mi) + vertex_steps[mi] * VN.row(mi);
      } else {  // new verts
        rV2mV[i] = finalV.size();
        finalV.push_back(rV.row(i));
        birthV.push_back(v);
      }
    }
    for (int i = 0; i < rF.rows(); i++) {
      finalF.emplace_back(
          Vec3i{rV2mV(rF(i, 0)), rV2mV(rF(i, 1)), rV2mV(rF(i, 2))});
      finalT.emplace_back(
          Vec4i{v, rV2mV(rF(i, 0)), rV2mV(rF(i, 1)), rV2mV(rF(i, 2))});
    }
    faces_to_erase.insert(VF[v].begin(), VF[v].end());
  }

  for (auto i : faces_to_erase) {
    finalF.erase(finalF.begin() + i);
  }
  {
    RowMati fT;
    RowMatd fV;
    RowMati fF;
    vec2eigen(finalV, fV);
    vec2eigen(finalF, fF);
    vec2eigen(finalT, fT);
    H5Easy::File file("complex_original_temp.h5", H5Easy::File::Overwrite);
    H5Easy::dump(file, "mV", mV);
    H5Easy::dump(file, "mF", mF);
    H5Easy::dump(file, "fV", fV);
    H5Easy::dump(file, "fF", fF);
    H5Easy::dump(file, "fT", fT);
    H5Easy::dump(file, "birth", birthV);
  }
}