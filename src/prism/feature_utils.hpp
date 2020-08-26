#pragma once

#include "common.hpp"
#include <list>
#include <set>
#include <vector>
namespace prism {
// reader from [Gao et al]
bool read_feature_graph(std::string path, Eigen::VectorXi &feat_nodes,
                        RowMati &feat_edges);

// split list of feature segments into disjoint chains.
// for cyclic loops, front == end;
bool split_feature_graph(const Eigen::VectorXi &feat_nodes,
                         const RowMati &feat_edges, int vnum,
                         std::vector<std::list<int>> &all_chain);

// segments the region adjacent to all chain verts.
// and split the regions into left/right, repectively for each chain.
void feature_chain_region_segments(const RowMatd &V, const RowMati &F,
                                   const std::vector<std::list<int>> chains,
                                   std::vector<std::set<int>> &feature_ignore);
} // namespace prism