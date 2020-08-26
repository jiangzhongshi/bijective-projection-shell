#ifndef PRISM_LOCAL_OPERATIONS_REMESH_PASS_HPP
#define PRISM_LOCAL_OPERATIONS_REMESH_PASS_HPP

#include "../common.hpp"
#include <list>
#include <vector>
#include <map>

class PrismCage;
namespace prism::local {
struct RemeshOptions {
  double distortion_bound = 0.1;
  double target_thickness = 0.01;
  bool parallel = true;
  double collapse_quality_threshold = 30;
  bool split_improve_quality = true;
  bool volume_centric = false; // volume quality etc.
  bool dynamic_hashgrid = false; // use a dynamic spatial hashgrid instead of static AABB

  std::function<double(const Vec3d &)> sizing_field;
  std::vector<double> target_adjustment;
  RemeshOptions() = default;
  RemeshOptions(int v_num, double edge_len) {
    sizing_field = [edge_len](const Vec3d &) { return edge_len; };
    target_adjustment.resize(v_num, 1);
  }
};
} // namespace prism::local

namespace prism::local {

constexpr auto shift_left = [](const auto &new_fid, const auto &new_shifts,
                               auto &F, auto &FF, auto &FFi) {
  constexpr auto roll_shift_left = [](auto &vec, int s) -> Vec3i {
    auto [a, b, c] = vec;
    if (s == 0)
      return {a, b, c};
    if (s == 1)
      return {b, c, a};
    else
      return {c, a, b};
  };
  for (int i = 0; i < new_shifts.size(); i++) {
    auto f = new_fid[i];
    auto s = new_shifts[i];
    // shift F,FF,FFi
    F[f] = roll_shift_left(F[f], s);
    FF[f] = roll_shift_left(FF[f], s);
    FFi[f] = roll_shift_left(FFi[f], s);
    // take care of the FFi for neighbors
    for (int j : {0, 1, 2}) {
      int f1 = FF[f][j];
      if (f1 == -1)
        continue;
      int e1 = FFi[f][j];
      FFi[f1][e1] = j;
    }
  }
};

int wildcollapse_pass(PrismCage &pc, RemeshOptions &);
void wildflip_pass(PrismCage &pc, const RemeshOptions &);
void wildsplit_pass(PrismCage &pc, RemeshOptions &);
void localsmooth_pass(PrismCage &pc, const RemeshOptions &);
void shellsmooth_pass(PrismCage &pc, const RemeshOptions &option);
} // namespace prism::local
#endif