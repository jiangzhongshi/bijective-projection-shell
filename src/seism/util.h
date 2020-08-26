#ifndef prism_UTIL
#define prism_UTIL
#include <Eigen/Core>
#include <utility>
#include <vector>
#include <list>

namespace prism::util{
template <typename DerivedF, typename Func, typename... Args>
inline int bfs_search_face(DerivedF&& F, DerivedF&& TT, Func succeed, size_t f0, Args&& ...args_for_succeed) {
  std::vector<bool> visited(TT.rows(), false);
  std::list<size_t> queue{f0};
  while (!queue.empty()) {
    // visit
    auto f = queue.front();
    queue.pop_front();
    visited[f] = true;

    if (succeed(F, TT, f))
        return f;

    // push neighbors
    for (int i = 0; i < 3; i++) {
      auto n = TT(f, i);
      if (n != -1 && !visited[n]) queue.push_back(n);
    }
  }
  return -1;
}
}

void split_boundary(
    const std::vector<int>& n,
    const std::vector<std::pair<int,int>>& P,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
);

std::tuple<double, double, bool> flip_avoiding_line_search(
  const Eigen::MatrixXi F,
  Eigen::MatrixXd& cur_v,
  Eigen::MatrixXd& dst_v,
  int max_iter,
  std::function<double(Eigen::MatrixXd&)> energy,
  double cur_energy);

#endif
