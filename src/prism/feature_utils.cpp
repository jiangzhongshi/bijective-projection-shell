#include "feature_utils.hpp"
#include "common.hpp"
#include <fstream>
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/spdlog.h>

#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>

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
  if (f.fail())
    return false;
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

bool prism::split_feature_graph(const Eigen::VectorXi &feat_nodes,
                                const RowMati &feat_edges, int vnum,
                                std::vector<std::list<int>> &all_chain) {
  RowMati feature_edges = feat_edges;
  // find corners
  std::map<int, int> valency;
  std::for_each(feature_edges.data(),
                feature_edges.data() + feature_edges.size(),
                [&valency](const auto &a) {
                  if (valency.find(a) == valency.end())
                    valency[a] = 0;
                  valency[a]++;
                });
  std::set<int> corners;
  for (auto [v, k] : valency)
    if (k != 2)
      corners.insert(v);
  std::for_each(feat_nodes.data(), feat_nodes.data() + feat_nodes.size(),
                [&corners](auto a) { corners.insert(a); });
  // split corners and build connectivity
  std::map<int, std::vector<int>> curve_conn;
  std::vector<int> extra_corner_ind;
  for (int i = 0, extra_vnum = 0; i < feature_edges.rows(); i++) {
    for (int j : {0, 1}) {
      auto &v = feature_edges(i, j);
      if (corners.find(v) !=
          corners.end()) { // replace each occurence with a fake one.
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
    if (k.size() == 1)
      assert(c >= vnum);
  }
#endif

  auto trace_through = [&curve_conn](int v0, int v1, auto f) {
    while (true) {
      auto itv1 = curve_conn.find(v1);
      if (itv1 == curve_conn.end())
        return true; // cycle
      auto e1 = itv1->second;
      assert(e1.size() > 0);
      curve_conn.erase(itv1);
      if (e1.size() < 2) // reached a corner.
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
    if (nb.size() == 2)
      chain.push_back(nb[1]);
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
                    if (a >= vnum)
                      a = extra_corner_ind[a - vnum];
                  });
    if (cycle) {
      chain.pop_back();
      assert(chain.front() == chain.back());
    }
    all_chain.emplace_back(chain);
  }
  return true;
}


void prism::feature_chain_region_segments(
    const RowMatd &V, const RowMati &F,
    const std::vector<std::list<int>> chains,
    std::vector<std::set<int>  /*refs*/> // 0 is left and 1 is right
       & feature_ignore) {
  std::set<int> feature_corners;
  for (auto &c : chains) {
    auto cf = c.front(), cb = c.back();
    if (cf == cb)
      continue;
    feature_corners.insert(cf);
    feature_corners.insert(cb);
  }
  Eigen::VectorXi VF, NI;
  RowMati TT, TTi;
  igl::vertex_triangle_adjacency(F, V.rows(), VF, NI);
  igl::triangle_triangle_adjacency(F, TT, TTi);
  for (int ci = 0; ci < chains.size(); ci++) {
    auto &c = chains[ci];
    std::set<int> faces_of_interest;
    std::set<int> seg0;
    std::map<std::pair<int, int>, int> edges_to_chain;
    for (auto it = std::next(c.begin()); it != c.end(); it++)
      edges_to_chain.emplace(std::pair(*std::prev(it), *it), ci);
    for (auto i : chains[ci]) {
      if (feature_corners.find(i) != feature_corners.end())
        continue; // corners are not added
      for (auto j = NI[i]; j < NI[i + 1]; j++)
        faces_of_interest.insert(VF[j]);
    }
    int flag = 0;
    std::function<void(int)> bfs_segment;
    bfs_segment = [&faces_of_interest, &seg0, &TT, &F, &flag, &bfs_segment, &ci,
                   &edges_to_chain](int fi) -> void {
      auto it = faces_of_interest.find(fi);
      if (it == faces_of_interest.end())
        return;
      faces_of_interest.erase(it);
      seg0.insert(*it);
      auto f = F.row(fi);
      for (auto j : {0, 1, 2}) {
        auto it0 = edges_to_chain.end();
        if (flag != -1)
          it0 = edges_to_chain.find({f[j], f[(j + 1) % 3]});
        if (it0 != edges_to_chain.end() && it0->second == ci) {
          assert(flag != -1);
          flag = 1; // the found it has a following edge (is a face on the
                    // left of chain)
          continue;
        }
        if (flag != 1)
          it0 = edges_to_chain.find({f[(j + 1) % 3], f[j]});
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
    if (flag == -1)
      std::swap(faces_of_interest, seg0);
    feature_ignore.emplace_back(seg0);              // left seg
    feature_ignore.emplace_back(faces_of_interest); // right seg
  }
  assert(feature_ignore.size() == 2*chains.size());
}