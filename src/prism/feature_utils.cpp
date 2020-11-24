#include "feature_utils.hpp"

#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/ostr.h>
#include "local_operations/mesh_coloring.hpp"
#include <spdlog/spdlog.h>

#include <fstream>

#include "common.hpp"

void prism::mark_feature_edges(const RowMatd &V, const RowMati &F,
                               double threshold, RowMati &feat) {
  RowMati TT, TTi;
  RowMatd FN;
  std::vector<std::array<int, 2>> feat_vec;
  igl::triangle_triangle_adjacency(F, TT, TTi);
  igl::per_face_normals(V, F, FN);
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      if (FN.row(i).dot(FN.row(TT(i, j))) < threshold) {
        auto v0 = F(i, j), v1 = F(i, (j + 1) % 3);
        feat_vec.push_back({std::min(v0, v1), std::max(v0, v1)});
      }
    }
  }
  std::sort(feat_vec.begin(), feat_vec.end());
  feat_vec.erase(std::unique(feat_vec.begin(), feat_vec.end()), feat_vec.end());
  feat = Eigen::Map<RowMati>(feat_vec[0].data(), feat_vec.size(), 2);
}

bool prism::read_feature_graph(std::string path, Eigen::VectorXi &feat_nodes,
                               RowMati &feat_edges) {
  // https://github.com/gaoxifeng/Feature-Preserving-Octree-Hex-Meshing/blob/67e6d59affb13116a260e40153a81f3bd3ac7256/io.cpp#L412
  // differently from obj, this is 0 based index, directly usable in the
  // program.
  struct Mesh_Feature {
    double angle_threshold;
    int orphan_curve;
    int orphan_curve_single;
    Eigen::VectorXi IN_corners;
    RowMati IN_v_pairs;
  };
  Mesh_Feature mf;
  std::fstream f(path, std::ios::in);
  if (f.fail()) return false;
  char s[1024];
  int cnum = 0, edgenum = 0;
  f.getline(s, 1023);
  sscanf(s, "%lf %i %i", &mf.angle_threshold, &mf.orphan_curve,
         &mf.orphan_curve_single);
  f.getline(s, 1023);
  sscanf(s, "%i %i", &cnum, &edgenum);
  mf.IN_corners.resize(cnum);
  for (int i = 0; i < cnum; i++) {
    f.getline(s, 1023);
    sscanf(s, "%i", &(mf.IN_corners[i]));
  }
  mf.IN_v_pairs.resize(edgenum, 2);
  for (int i = 0; i < edgenum; i++) {
    f.getline(s, 1023);
    int v0 = -1, v1 = -1;
    sscanf(s, "%i %i", &v0, &v1);
    mf.IN_v_pairs(i, 0) = v0;
    mf.IN_v_pairs(i, 1) = v1;
  }
  f.close();
  feat_nodes = std::move(mf.IN_corners);
  feat_edges = std::move(mf.IN_v_pairs);
  return true;
}

bool prism::feature_chains_from_edges(const Eigen::VectorXi &feat_nodes,
                                      const RowMati &feat_edges, int vnum,
                                      std::vector<std::list<int>> &all_chain) {
  RowMati feature_edges = feat_edges;
  // find corners
  std::map<int, int> valency;
  std::for_each(feature_edges.data(),
                feature_edges.data() + feature_edges.size(),
                [&valency](const auto &a) {
                  if (valency.find(a) == valency.end()) valency[a] = 0;
                  valency[a]++;
                });
  std::set<int> corners;
  for (auto [v, k] : valency)
    if (k != 2) corners.insert(v);
  std::for_each(feat_nodes.data(), feat_nodes.data() + feat_nodes.size(),
                [&corners](auto a) { corners.insert(a); });
  // split/duplicate corners and build connectivity
  // curve connectivity is defined by v0 -> list of adjacent (up to two).
  std::map<int, std::vector<int>> curve_conn;
  std::vector<int> extra_corner_ind;
  for (int i = 0, extra_vnum = 0; i < feature_edges.rows(); i++) {
    for (int j : {0, 1}) {
      auto &v = feature_edges(i, j);
      if (corners.find(v) !=
          corners.end()) {  // replace each occurence with a fake one.
        extra_corner_ind.push_back(v);
        v = vnum + extra_vnum++;
      }
      if (curve_conn.find(v) == curve_conn.end())
        curve_conn.emplace(v, std::vector<int>{});
    }
    auto v0 = feature_edges(i, 0);
    auto v1 = feature_edges(i, 1);
    assert(corners.find(v0) == corners.end());
    assert(corners.find(v1) == corners.end());
    curve_conn[v0].push_back(v1);
    curve_conn[v1].push_back(v0);
  }

#ifndef NDEBUG
  for (auto &[c, k] : curve_conn) {
    assert(k.size() <= 2);
    if (k.size() == 1) assert(c >= vnum);
  }
#endif

  auto trace_through = [&curve_conn](int v0, int v1, auto f) {
    while (true) {
      auto itv1 = curve_conn.find(v1);
      if (itv1 == curve_conn.end()) return true;  // cycle
      auto e1 = itv1->second;
      assert(e1.size() > 0);
      curve_conn.erase(itv1);
      if (e1.size() < 2)  // reached a corner.
        return false;
      auto tmp = v1;
      assert(e1[0] == v0 || e1[1] == v0);
      v1 = e1[0] + e1[1] - v0;
      v0 = tmp;
      *f = (v1);
    }
  };

  while (curve_conn.size() > 0) {
    auto it = curve_conn.begin();
    auto nb = it->second;
    auto v0 = it->first, v1 = nb[0];
    std::list<int> chain{v1, v0};
    if (nb.size() == 2) chain.push_back(nb[1]);
    curve_conn.erase(it);
    bool cycle = trace_through(v0, v1, std::front_inserter(chain));
    if (!cycle && nb.size() == 2) {
      assert(chain.front() >= vnum);
      v1 = nb[1];
      bool c = trace_through(v0, v1, std::back_inserter(chain));
      assert(!c);
      assert(chain.back() >= vnum);
    }
    std::for_each(chain.begin(), chain.end(),
                  [&extra_corner_ind, vnum](auto &a) {
                    if (a >= vnum) a = extra_corner_ind[a - vnum];
                  });
    if (cycle) {
      chain.pop_back();
      assert(chain.front() == chain.back());
    }
    all_chain.emplace_back(chain);
  }
  return true;
}

void feature_chain_region_segments_legacy(
    const RowMati &F, int vnum, const std::vector<std::list<int>> chains,
    std::vector<std::set<int> /*refs*/>  // 0 is left and 1 is right
        &feature_side_region) {
  std::set<int> feature_corners;
  for (auto &c : chains) {
    auto cf = c.front(), cb = c.back();
    if (cf == cb) continue;
    feature_corners.insert(cf);
    feature_corners.insert(cb);
  }
  Eigen::VectorXi VF, NI;
  RowMati TT, TTi;
  igl::vertex_triangle_adjacency(F, vnum, VF, NI);
  igl::triangle_triangle_adjacency(F, TT, TTi);
  for (int ci = 0; ci < chains.size(); ci++) {
    auto &c = chains[ci];
    std::set<int> faces_of_interest;
    std::set<int> seg0;
    std::map<std::pair<int, int>, int> edges_to_chain;
    for (auto it = std::next(c.begin()); it != c.end(); it++)
      edges_to_chain.emplace(std::pair(*std::prev(it), *it), ci);
    // first, collect all adjacent faces to the current chain
    for (auto i : chains[ci]) {
      if (feature_corners.find(i) != feature_corners.end())
        continue;  // corners are not added
      for (auto j = NI[i]; j < NI[i + 1]; j++) faces_of_interest.insert(VF[j]);
    }

    // bfs coloring, flag is for swapping left/right.
    int flag = 0;
    std::function<void(int)> bfs_segment;
    bfs_segment = [&faces_of_interest, &seg0, &TT, &F, &flag, &bfs_segment, &ci,
                   &edges_to_chain](int fi) -> void {
      auto it = faces_of_interest.find(fi);
      if (it == faces_of_interest.end()) return;
      seg0.insert(*it);
      faces_of_interest.erase(it);
      auto f = F.row(fi);
      for (auto j : {0, 1, 2}) {
        auto it0 = edges_to_chain.end();
        if (flag != -1) it0 = edges_to_chain.find({f[j], f[(j + 1) % 3]});
        if (it0 != edges_to_chain.end() && it0->second == ci) {
          assert(flag != -1);
          flag = 1;  // the found it has a following edge (is a face on the
                     // left of chain)
          continue;
        }
        if (flag != 1) it0 = edges_to_chain.find({f[(j + 1) % 3], f[j]});
        if (it0 != edges_to_chain.end() && it0->second == ci) {
          assert(flag != 1);
          flag = -1;
          continue;
        }
        bfs_segment(TT(fi, j));
      }
    };
    bfs_segment(*faces_of_interest.begin());
    assert(flag != 0);
    if (flag == -1) std::swap(faces_of_interest, seg0);
    feature_side_region.emplace_back(seg0);               // left seg
    feature_side_region.emplace_back(faces_of_interest);  // right seg
  }
  assert(feature_side_region.size() == 2 * chains.size());
}

void prism::feature_chain_region_segments(
    const RowMati &F, int vnum, const std::vector<std::list<int>> chains,
    std::vector<std::set<int> /*refs*/>  // 0 is left and 1 is right
        &feature_side_region,
        std::vector<std::set<int>>& region_around_chains) {
  std::vector<std::vector<int>> VF, VFi;
  igl::vertex_triangle_adjacency(vnum, F, VF, VFi);
  RowMati TT, TTi;
  igl::triangle_triangle_adjacency(F, TT, TTi);
  for (auto ch : chains) {
    for (auto it = ch.begin(); std::next(it) != ch.end(); it++) {
      auto v0 = *it, v1 = *std::next(it);
      auto &nb = VF[v0], &nbi = VFi[v0];
      for (int i = 0; i < nb.size(); i++) {
        auto f = nb[i];
        assert(F(f, nbi[i]) == v0);
        if (F(f, (nbi[i] + 1) % 3) == v1) {  // lefty
          auto &f1 = TT(f, nbi[i]), &e1 = TTi(f, nbi[i]);
          TT(f1, e1) = -1;
          TTi(f1, e1) = -1;
          f1 = -1;
          e1 = -1;
        }
      }
    }
  }
  for (int ci = 0; ci < chains.size(); ci++) {
    auto &ch = chains[ci];
    std::set<int> faces_of_interest;
    std::set<int> seg0;
    std::map<std::pair<int, int>, int> edges_to_chain;
    for (auto it = std::next(ch.begin()); it != ch.end(); it++)
      edges_to_chain.emplace(std::pair(*std::prev(it), *it), ci);
    // first, collect all adjacent faces to the current chain
    for (auto i : chains[ci]) {
      faces_of_interest.insert(VF[i].begin(), VF[i].end());
    }

    // bfs coloring
    // initial coloring.
    std::set<int> solid_left, solid_right;
    for (auto it = ch.begin(); std::next(it) != ch.end(); it++) {
      auto v0 = *it, v1 = *std::next(it);
      auto &nb = VF[v0], &nbi = VFi[v0];
      for (int i = 0; i < nb.size(); i++) {
        auto f = nb[i];
        assert(F(f, nbi[i]) == v0);
        if (F(f, (nbi[i] + 1) % 3) == v1) {  // lefty
          solid_left.insert(f);
        }
        if (F(f, (nbi[i] + 2) % 3) == v1) {  // righty
          solid_right.insert(f);
        }
      }
    }

    std::function<void(int, const std::set<int> &)> bfs_segment;
    std::set<int> remain_faces;
    bfs_segment = [&remain_faces, &TT, &F, &bfs_segment](
                      int fi, const std::set<int> &solid_reject) -> void {
      if (fi < 0 || solid_reject.find(fi) != solid_reject.end()) return;
      auto it = remain_faces.find(fi);
      if (it == remain_faces.end()) {
        return;
      }
      remain_faces.erase(it);
      for (auto j : {0, 1, 2}) {
        bfs_segment(TT(fi, j), solid_reject);
      }
    };
    remain_faces = faces_of_interest;
    bfs_segment(*solid_left.begin(), solid_right);
    feature_side_region.emplace_back(remain_faces);  // rejector for the left
    remain_faces = faces_of_interest;
    bfs_segment(*solid_right.begin(), solid_left);
    feature_side_region.emplace_back(remain_faces);  // rejector for the right
    region_around_chains.emplace_back(faces_of_interest);
  }
  assert(feature_side_region.size() == 2 * chains.size());
}

std::vector<std::list<int>> prism::recover_chains_from_meta_edges(
    const std::map<std::pair<int, int>, std::pair<int, std::vector<int>>>
        &meta) {
  auto glue_together = [](std::vector<std::vector<int>> some_segs) {
    std::map<int, std::list<int>> v0_seg;
    for (auto s : some_segs) {
      v0_seg.emplace(s[0], std::list<int>(s.begin(), s.end()));
    }
    auto it = v0_seg.begin();
    while (v0_seg.size() > 1) {
      auto v1 = it->second.back(), v0 = it->first;
      auto it1 = v0_seg.find(v1);
      if (it1 != v0_seg.end()) {
        assert(it1->second.size() >= 2);
        it->second.insert(it->second.end(), std::next(it1->second.begin()),
                          it1->second.end());
        v0_seg.erase(it1);
      } else {
        it++;
      }
    }
    return v0_seg.begin()->second;
  };
  std::vector<std::vector<std::vector<int>>> collect;
  for (auto [m, cid_chain] : meta) {
    auto [cid, chain] = cid_chain;
    auto n = chain.size();
    if (collect.size() <= cid) collect.resize(cid + 1);
    collect[cid].push_back(chain);
  }
  std::vector<std::list<int>> chains(collect.size());
  for (int i = 0; i < collect.size(); i++) {
    chains[i] = glue_together(collect[i]);
  }
  return chains;
}

auto vv2fe = [](auto &v0, auto &v1, const auto &F, const auto &VF) {
  std::vector<int> common;
  std::set_intersection(VF[v0].begin(), VF[v0].end(), VF[v1].begin(),
                        VF[v1].end(), std::back_inserter(common));
  for (auto f : common) {
    for (int j = 0; j < 3; j++) {
      if (F(f, j) == v0 && F(f, (j + 1) % 3) == v1) {
        return std::pair(f, j);
      }
    }
  }
  return std::pair(-1, -1);
};

std::tuple<RowMatd, RowMati, RowMati> prism::subdivide_feature_triangles(const RowMatd &mV, const RowMati &mF,
                                 const RowMati &feature_edges) {
  // splits
  RowMati FF, FFi;
  igl::triangle_triangle_adjacency(mF, FF, FFi);
  Eigen::VectorXi colors = Eigen::VectorXi::Zero(mF.rows());
  std::vector<std::vector<int>> VF, VFi;
  igl::vertex_triangle_adjacency(mV.rows(), mF, VF, VFi);
  // initialize red (2)
  RowMati feature_fe(mF.rows(), 3);
  feature_fe.setZero();
  for (auto e = 0; e < feature_edges.rows(); e++) {
    auto [f, j] = vv2fe(feature_edges(e, 0), feature_edges(e, 1), mF, VF);
    feature_fe(f, j) = 1;
    colors[f] = 2;
    colors[FF(f, j)] = 2;
  }
  RowMati edge_vert = -RowMati::Ones(mF.rows(), 3);
  prism::local::red_green_coloring(mF, FF, colors);
  spdlog::set_level(spdlog::level::debug);
  std::vector<Vec3d> V;
  std::vector<Vec3i> F;
  eigen2vec(mV, V);
  eigen2vec(mF, F);
  for (int i = 0; i < mF.rows(); i++) {
    if (colors[i] != 1) continue;  // only deal with green here
    int e = [&colors, &FF, &i]() {
      for (int j = 0; j < 3; j++) {
        if (FF(i, j) != -1 && colors[FF(i, j)] == 2) {
          return j;
        }  // find red neighbor
      }
      return -1;
    }();
    assert(e != -1);

    int ux = V.size();
    V.push_back((V[F[i][e]] + V[F[i][(e + 1) % 3]]) / 2);
    edge_vert(i, e) = ux;
    edge_vert(FF(i, e), FFi(i, e)) = ux;
    auto v0 = mF(i, e), v1 = mF(i, (e + 1) % 3), v2 = mF(i, (e + 2) % 3);
    F.emplace_back(Vec3i{v0, ux, v2});
    F.emplace_back(Vec3i{v1, v2, ux});
  }

  std::vector<std::array<int, 2>> sub_features;
  for (int i = 0; i < colors.size(); i++) {
    if (colors[i] != 2) continue;
    for (int e = 0; e < 3; e++) {
      if (edge_vert(i, e) == -1) {  // assigned
        edge_vert(i, e) = edge_vert(FF(i, e), FFi(i, e)) = V.size();
        V.push_back((V[mF(i, e)] + V[mF(i, (e + 1) % 3)]) / 2);
      }
      if (feature_fe(i, e)) {
        sub_features.push_back({mF(i, e), edge_vert(i, e)});
        sub_features.push_back({edge_vert(i, e), mF(i, (e + 1) % 3)});
      }
    }
    auto v0 = mF(i, 0), v1 = mF(i, 1), v2 = mF(i, 2);
    auto e0 = edge_vert(i, 0), e1 = edge_vert(i, 1), e2 = edge_vert(i, 2);
    F.emplace_back(Vec3i{v0, e0, e2});
    F.emplace_back(Vec3i{v1, e1, e0});
    F.emplace_back(Vec3i{v2, e2, e1});
    F.emplace_back(Vec3i{e0, e1, e2});
  }  // insert red faces
  for (int i = 0; i < colors.size(); i++)
    if (colors[i] > 0) {  // removed
      F[i] = F.back();
      F.pop_back();
    }

  return std::tuple(RowMatd(Eigen::Map<RowMatd>(V[0].data(), V.size(), 3)),
                    RowMati(Eigen::Map<RowMati>(F[0].data(), F.size(), 3)),
                    RowMati(Eigen::Map<RowMati>(sub_features[0].data(),
                                                sub_features.size(), 2)));
}