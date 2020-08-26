#ifndef PRISM_SPATIAL_HASH_AABB_HASH
#define PRISM_SPATIAL_HASH_AABB_HASH
#include <list>
#include <memory>
#include <set>
/// @brief An entry into the hash grid as a (key, value) pair.

#include "../common.hpp"
namespace GEO {
class Box;
};
namespace prism {
struct HashItem {
  int key; /// @brief The key of the item.
  int id;  /// @brief The value of the item.
  std::shared_ptr<GEO::Box>
      aabb; /// @brief The axis-aligned bounding box of the element
  HashItem(int key, int id) : key(key), id(id) {}

  /// @brief Compare HashItems by their keys for sorting.
  bool operator<(const HashItem &other) const;
};

using HashMap =
    std::map<std::tuple<int, int, int>, std::shared_ptr<std::list<int>>>;
using HashPtr =
    std::pair<std::shared_ptr<std::list<int>>, std::list<int>::iterator>;
struct HashGrid {
  HashGrid(const Vec3d &lower, const Vec3d &upper, double cell)
      : m_domain_min(lower), m_domain_max(upper), m_cell_size(cell) {
    m_grid_size = int(std::ceil((upper - lower).maxCoeff() / m_cell_size));
  };
  HashGrid(const std::vector<Vec3d> &V, const std::vector<Vec3i> &F);
  void insert_triangles(const std::vector<Vec3d> &V,
                        const std::vector<Vec3i> &F,
                        const std::vector<int> &fid);

  bool self_intersect() const;
  void query(const Vec3d &lower, const Vec3d &upper, std::set<int> &) const;
  void add_element(const Vec3d &lower, const Vec3d &upper, const int index);
  void remove_element(const int index);
  void bound_convert(const Vec3d &, Eigen::Array3i &) const;

  // spatial partition parameters
  Vec3d m_domain_min;
  Vec3d m_domain_max;
  double m_cell_size;
  size_t m_grid_size;

  HashMap m_face_items;
  std::vector<std::vector<HashPtr>> face_stores; // facilitate element removal
};
} // namespace prism
#endif