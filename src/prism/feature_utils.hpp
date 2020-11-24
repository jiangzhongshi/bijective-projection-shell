#ifndef PRISM_FEATUE_UTILS_HPP
#define PRISM_FEATUE_UTILS_HPP

#include <list>
#include <set>
#include <vector>

#include "common.hpp"
namespace prism {
// reader from [Gao et al]
bool read_feature_graph(std::string path, Eigen::VectorXi &feat_nodes,
                        RowMati &feat_edges);

// compute the dot product of adjacent face normals (TT and FN)
// and mark v0-v1 if below threshold.
void mark_feature_edges(const RowMatd &V, const RowMati &F, double thre,
                        RowMati &feat);
// connect list of feature edges, and split into disjoint chains.
// feat_nodes is optional and high valence connectors are automatically
// detected. Note: Cyclic loops is not specially marked, since front == end;
bool feature_chains_from_edges(const Eigen::VectorXi &feat_nodes,
                               const RowMati &feat_edges, int vnum,
                               std::vector<std::list<int>> &all_chain);

// mark the region around the chains
// for each chain, mark the reference that should be reject by left/right.
void feature_chain_region_segments(
    const RowMati &F, int vnum, const std::vector<std::list<int>> chains,
    std::vector<std::set<int>> &feature_ignore,
    std::vector<std::set<int>> &region_around_chain);

std::vector<std::list<int>> recover_chains_from_meta_edges(
    const std::map<std::pair<int, int>, std::pair<int, std::vector<int>>>
        &meta);

// a helper function for testing feature snapping.
// given input triangle mesh with feature annotation, subdivde (red-green) the
// triangles along the feature creating more challenge for curve only snapping.
std::tuple<RowMatd, RowMati, RowMati> subdivide_feature_triangles(
    const RowMatd &mV, const RowMati &mF, const RowMati &feature_edges);
}  // namespace prism

#endif