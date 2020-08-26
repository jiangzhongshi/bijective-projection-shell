#include "AABB_hash.hpp"
#include <geogram/basic/geometry.h>
#include <igl/avg_edge_length.h>
#include <numeric>

bool prism::HashItem::operator<(const prism::HashItem &other) const {
  return std::tie(key, id) < std::tie(other.key, other.id);
}

prism::HashGrid::HashGrid(const std::vector<Vec3d> &V,
                          const std::vector<Vec3i> &F) {
  RowMatd matV;
  RowMati matF;
  vec2eigen(V, matV);
  vec2eigen(F, matF);
  m_domain_min = matV.colwise().minCoeff();
  m_domain_max = matV.colwise().maxCoeff();
  double avg_len = igl::avg_edge_length(matV, matF);
  m_cell_size = 2 * avg_len;
  m_grid_size =
      int(std::ceil((m_domain_max - m_domain_min).maxCoeff() / m_cell_size));
  std::vector<int>fid(F.size());
  std::iota(fid.begin(), fid.end(), 0);
  insert_triangles(V, F, fid);
}

void prism::HashGrid::bound_convert(const Vec3d &aabb_min,
                                    Eigen::Array3i &int_min) const {
  int_min = ((aabb_min - m_domain_min) / m_cell_size).cast<int>().array();
  assert((int_min.minCoeff() >= -1));
  assert((int_min.maxCoeff() <= m_grid_size));
  int_min.max(0).min(m_grid_size - 1);
}

void prism::HashGrid::add_element(const Vec3d &aabb_min, const Vec3d &aabb_max,
                                  const int index) {
  assert(face_stores[index].size() == 0);
  Eigen::Array3i int_min, int_max;
  bound_convert(aabb_min, int_min);
  bound_convert(aabb_max, int_max);
  assert(int_min[0] <= int_max[0]);

  auto &hg = m_face_items;
  for (int x = int_min.x(); x <= int_max.x(); ++x)
    for (int y = int_min.y(); y <= int_max.y(); ++y)
      for (int z = int_min.z(); z <= int_max.z(); ++z) {
        auto key = std::tuple(x, y, z);
        auto it = hg.lower_bound(key);
        if (it != hg.end() && it->first == key) {
          it->second->push_back(index);
        } else {
          it = hg.emplace_hint(
              it, key, std::make_shared<std::list<int>>(std::list<int>{index}));
        }
        auto ptr = std::prev(it->second->end());
        face_stores[index].emplace_back(it->second, ptr);
      }
}

void prism::HashGrid::remove_element(const int index) {
  for (auto [ptr_list, iter_list] : face_stores[index]) {
    ptr_list->erase(iter_list);
  }
  face_stores[index].clear();
}

void prism::HashGrid::query(const Vec3d &aabb_min, const Vec3d &aabb_max,
                            std::set<int> &result) const {
  Eigen::Array3i int_min, int_max;
  bound_convert(aabb_min, int_min);
  bound_convert(aabb_max, int_max);
  assert(int_min[0] <= int_max[0]);

  auto &hg = m_face_items;
  std::vector<int> vec_result;
  for (int x = int_min.x(); x <= int_max.x(); ++x)
    for (int y = int_min.y(); y <= int_max.y(); ++y)
      for (int z = int_min.z(); z <= int_max.z(); ++z) {
        auto key = std::tuple(x, y, z);
        auto it = hg.find(key);
        if (it != hg.end()) {
          vec_result.insert(vec_result.end(), it->second->begin(),
                            it->second->end());
        }
      }
  result = std::set<int>(vec_result.begin(), vec_result.end());
  // std::set(query_result.begin(), query_result.end());
  // std::unique(query_result.begin(), query_result.end());
}

void prism::HashGrid::insert_triangles(const std::vector<Vec3d> &V,
                                       const std::vector<Vec3i> &F,
                                       const std::vector<int> &fid) {
  face_stores.resize(F.size());
  for (auto i : fid) {
    auto [v0, v1, v2] = F[i];
    Eigen::Matrix3d local;
    for (auto k : {0, 1, 2})
      local.row(k) = V[F[i][k]];
    auto aabb_min = local.colwise().minCoeff();
    auto aabb_max = local.colwise().maxCoeff();
    add_element(aabb_min, aabb_max, i);
  }

#ifndef NDEBUG
  for (auto i : fid) {
    for (auto [ptr_list, iter_list] : face_stores[i]) {
      assert(*iter_list == i);
    }
  }
#endif
}

bool prism::HashGrid::self_intersect() const{
  for (auto [_, l_ptr]: m_face_items) {
    auto n = l_ptr->size();
    for (int i=0; i<n;i++)
      for(int j=i+1; j<n;j++) { // test i-j pair

      }
  }
  assert(false);
}