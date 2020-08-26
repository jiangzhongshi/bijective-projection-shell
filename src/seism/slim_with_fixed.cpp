#include <igl/slim.h>
#include "slim_with_fixed.h"
#include <igl/min_quad_with_fixed.h>

bool solve_weighted_arap_with_hard(SLIMFixedData& s,
                                    const Eigen::MatrixXd &V,
                                    const Eigen::MatrixXi &F,
                                    Eigen::MatrixXd &uv,
                                    Eigen::VectorXi &soft_b_p,
                                    Eigen::MatrixXd &soft_bc_p)
{
  using namespace Eigen;

  Eigen::SparseMatrix<double> L;
  igl::slim::build_linear_system(s,L);

  // Eigen::SparseMatrix<double> LmLt = (L-L.transpose());
  // std::cout<<"InSymmetric"<<(Eigen::Map<const Eigen::VectorXd>(LmLt.valuePtr(), LmLt.size()).maxCoeff())<<std::endl;
  if (!igl::is_symmetric(L, 1.0))
   {
     std::cout<<"SLIM warning: returning because L is not Symmetric"<<std::endl;
    return false;
   }
    // L = (L + Eigen::SparseMatrix<double>(L.transpose()))/2;

  igl::Timer t;
  
  Eigen::VectorXd Uc;
  assert(s.dim == 2);
#ifndef CHOLMOD
  Eigen::VectorXi known(s.hard_constraints.rows()*2);
  
  for (int i=0; i<s.hard_constraints.rows(); i++) {
    known(i) = s.hard_constraints(i);
    known(i + s.hard_constraints.rows()) = s.hard_constraints(i) + s.v_num;
  }
  VectorXd flat_uv  = Eigen::Map<Eigen::VectorXd>(s.V_o.data(), s.v_num*2);

  Eigen::VectorXd Y = igl::slice(flat_uv, known);

  Eigen::MatrixXd B = -s.rhs;

  Eigen::SparseMatrix<double> Aeq;
  Eigen::MatrixXd Beq;
  bool status = igl::min_quad_with_fixed(L, B, known, Y, Aeq, Beq, true, Uc);
  if (!status)
   {
     std::cout<<"min_quad_with_fixed failed"<<std::endl;
    return false;
   }
#else
    CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    Uc = solver.compute(L).solve(s.rhs);
#endif
  for (int i = 0; i < s.dim; i++)
    uv.col(i) = Uc.block(i * s.v_n, 0, s.v_n, 1);
  return true;

}


Eigen::MatrixXd slim_fixed_solve(SLIMFixedData &data, int iter_num)
{
  for (int i = 0; i < iter_num; i++)
  {
    Eigen::MatrixXd dest_res;
    dest_res = data.V_o;

    // Solve Weighted Proxy
    igl::slim::update_weights_and_closest_rotations(data, dest_res);
    if (!solve_weighted_arap_with_hard(data,data.V, data.F, dest_res, data.b, data.bc))
      return data.V_o;

    double old_energy = data.energy;

    std::function<double(Eigen::MatrixXd &)> compute_energy = [&](
        Eigen::MatrixXd &aaa) { return igl::slim::compute_energy(data,aaa); };

    data.energy = igl::flip_avoiding_line_search(data.F, data.V_o, dest_res, compute_energy,
                                                 data.energy * data.mesh_area) / data.mesh_area;
  }
  return data.V_o;
}
