#ifndef PRISM_PRISMCAGE_HPP
#define PRISM_PRISMCAGE_HPP

#include <array>
#include <list>
#include <memory>
#include <set>
#include <tuple>

#include "common.hpp"

namespace prism::geogram {
class AABB;
};
namespace prism {
struct HashGrid;
};

struct PrismCage {
  enum class SeparateType { kShell, kSurface, kNone };
  ///////////////////////////////////////
  // The input reference surface, with fixed structure throughout the
  // optimization. contains V, F and acceleration AABB.
  ///////////////////////////////////////
  struct RefSurf {
    RowMatd snapV;
    RowMatd V;
    RowMati F;
    std::unique_ptr<prism::geogram::AABB> aabb;
  };
  RefSurf ref;

  void serialize(std::string filename);

  ///////////////////////////////////////
  // Data for the Cage
  // Base Vertex, Top Vertex, F, TT, TTi
  // All with std::vector to enable dynamic change.
  ///////////////////////////////////////
  std::vector<Vec3d> base;
  std::vector<Vec3d> top;
  std::vector<Vec3d> mid;
  std::vector<Vec3i> F;

  std::vector<std::set<int>> track_ref;
  std::shared_ptr<prism::HashGrid> base_grid = nullptr;
  std::shared_ptr<prism::HashGrid> top_grid = nullptr;

  Eigen::VectorXi feature_corners;
  std::vector<std::list<int>> feature_chains;
  RowMati feature_edges;
  std::map<std::pair<int, int>, std::pair<int, std::vector<int>>> meta_edges;
  ///////////////////////////////////////
  // Constructor and Initialize cage
  ///////////////////////////////////////
  PrismCage() = default;
  PrismCage(std::string);
  PrismCage(const RowMatd &vert, const RowMati &face, double dooseps = 0.2,
            double initial_step = 1e-4, SeparateType st=SeparateType::kSurface,
            const Eigen::VectorXi &feat = Eigen::VectorXi(),
            const RowMati &feat_vv = RowMati());
  void load_from_hdf5(std::string);
  void construct_cage(const RowMatd &);
  void init_track();
  void cleanup_empty_faces(Eigen::VectorXi &NI, Eigen::VectorXi &NJ);
};

#endif