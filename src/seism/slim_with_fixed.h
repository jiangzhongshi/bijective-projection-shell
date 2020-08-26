#ifndef SLIM_WITH_FIXED_H
#define SLIM_WITH_FIXED_H

struct SLIMFixedData: igl::SLIMData 
{
  Eigen::VectorXi hard_constraints;
};
Eigen::MatrixXd slim_fixed_solve(SLIMFixedData &data, int iter_num);
#endif