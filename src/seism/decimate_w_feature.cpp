#include <igl/decimate.h>
#include <igl/is_edge_manifold.h>
#include <igl/shortest_edge_and_midpoint.h>

using Xi = Eigen::MatrixXi;
using Vi = Eigen::VectorXi;
using Vd = Eigen::VectorXd;
using Xd = Eigen::MatrixXd;

bool decimate_with_feature(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                          int cons_num, Eigen::MatrixXd& U, Eigen::MatrixXi& G, Eigen::VectorXi& I, Eigen::VectorXi& J) {
  // Original number of faces
  const int orig_m = F.rows();
  // Tracking number of faces
  int m = F.rows();
  typedef Eigen::MatrixXd DerivedV;
  typedef Eigen::MatrixXi DerivedF;
  DerivedV VO;
  DerivedF FO;
  igl::connect_boundary_to_infinity(V,F,VO,FO);
  // decimate will not work correctly on non-edge-manifold meshes. By extension
  // this includes meshes with non-manifold vertices on the boundary since these
  // will create a non-manifold edge when connected to infinity.
  if(!igl::is_edge_manifold(FO))
  {
    return false;
  }

  auto constraint_stopping_condition = [&cons_num](const int e,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & /*F*/,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & /*EMAP*/,
  const Eigen::MatrixXi & /*EF*/,
  const Eigen::MatrixXi & /*EI*/,
  double & cost,
  Eigen::RowVectorXd & p)->void{
      int v0=E(e,0), v1 = E(e,1);
      if (v0<cons_num && v1 < cons_num){
        cost = std::numeric_limits<double>::infinity();
        p = V.row(v0);
        return ;
      }
    for(auto&v :{v0,v1})
        if (v<cons_num){
        cost = (V.row(E(e,0))-V.row(E(e,1))).norm();
        p = V.row(v);
      }
      cost = (V.row(E(e,0))-V.row(E(e,1))).norm();
      p = 0.5*(V.row(E(e,0))+V.row(E(e,1)));
  };
  bool ret = igl::decimate(
    VO, FO, constraint_stopping_condition,
    igl::max_faces_stopping_condition(m,orig_m,0),
    U, G, J, I);
  const Eigen::Array<bool,Eigen::Dynamic,1> keep = (J.array()<orig_m);
  igl::slice_mask(Eigen::MatrixXi(G),keep,1,G);
  igl::slice_mask(Eigen::VectorXi(J),keep,1,J);
  Eigen::VectorXi _1,I2;
  igl::remove_unreferenced(Eigen::MatrixXd(U),Eigen::MatrixXi(G),U,G,_1,I2);
  igl::slice(Eigen::VectorXi(I),I2,1,I);
  return ret;
}


#include <map>
#include <utility>

bool distance_preserve_decimate(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                          int cons_num, Eigen::MatrixXd& U, Eigen::MatrixXi& G, Eigen::VectorXi& I, Eigen::VectorXi& J,
                          std::map<std::pair<int,int>, double>& dists) {
  auto sorted_pair = [](int v0, int v1){return std::make_pair(std::min(v0,v1), std::max(v0,v1));};
  // Original number of faces
  const int orig_m = F.rows();
  // Tracking number of faces
  int m = F.rows();
  typedef Eigen::MatrixXd DerivedV;
  typedef Eigen::MatrixXi DerivedF;
  DerivedV VO;
  DerivedF FO;
  igl::connect_boundary_to_infinity(V,F,VO,FO);
  // decimate will not work correctly on non-edge-manifold meshes. By extension
  // this includes meshes with non-manifold vertices on the boundary since these
  // will create a non-manifold edge when connected to infinity.
  if(!igl::is_edge_manifold(FO))
  {
    return false;
  }
  
  for (int i=0; i<F.rows(); i++) {
    for (int j=0; j<3; j++) {
      int v0 = F(i,j);
      int v1 = F(i,(j+1)%3);
      double cost = (V.row(v0)-V.row(v1)).norm();
      auto vpair = std::make_pair(std::min(v0,v1), std::max(v0,v1));

      dists[vpair] = cost;
    }
  }
  int recent_collapsed_vert = -1;
  auto constraint_stopping_condition = [&cons_num, &dists, &sorted_pair, &recent_collapsed_vert](const int e,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & /*F*/,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & /*EMAP*/,
  const Eigen::MatrixXi & /*EF*/,
  const Eigen::MatrixXi & /*EI*/,
  double & cost,
  Eigen::RowVectorXd & p)->void{
      int v0=E(e,0), v1 = E(e,1);
      auto vpair = std::make_pair(std::min(v0,v1), std::max(v0,v1));
      if (dists.find(vpair) == dists.end()) { // not found, need to recompute
        auto d0 = dists.find(sorted_pair(v0,recent_collapsed_vert));
        auto d1 = dists.find(sorted_pair(v1,recent_collapsed_vert));

        dists[vpair] = (d0->second) + (d1->second);
        // dists.erase(d0);
        // dists.erase(d1);

        p = V.row(0); 
      }
      else{
        cost = dists[vpair];
        p = V.row(0);
      }
      if (v0<cons_num && v1 < cons_num) // never collapse two constraints together
        cost = std::numeric_limits<double>::infinity();
  };
  const auto always_try = [&recent_collapsed_vert](
    const Eigen::MatrixXd &                                         ,/*V*/
    const Eigen::MatrixXi &                                         ,/*F*/
    const Eigen::MatrixXi & E                                        ,/*E*/
    const Eigen::VectorXi &                                         ,/*EMAP*/
    const Eigen::MatrixXi &                                         ,/*EF*/
    const Eigen::MatrixXi &                                         ,/*EI*/
    const std::set<std::pair<double,int> > &                        ,/*Q*/
    const std::vector<std::set<std::pair<double,int> >::iterator > &,/*Qit*/
    const Eigen::MatrixXd &                                         ,/*C*/
    const int e                                                       /*e*/
    ) -> bool { 
      recent_collapsed_vert = std::max(E(e,0), E(e,1));
      return true;};

  const auto never_care = [](
  const Eigen::MatrixXd &                                         ,   /*V*/
  const Eigen::MatrixXi &                                         ,   /*F*/
  const Eigen::MatrixXi &                                         ,   /*E*/
  const Eigen::VectorXi &                                         ,/*EMAP*/
  const Eigen::MatrixXi &                                         ,  /*EF*/
  const Eigen::MatrixXi &                                         ,  /*EI*/
  const std::set<std::pair<double,int> > &                        ,   /*Q*/
  const std::vector<std::set<std::pair<double,int> >::iterator > &, /*Qit*/
  const Eigen::MatrixXd &                                         ,   /*C*/
  const int                                                       ,   /*e*/
  const int                                                       ,  /*e1*/
  const int                                                       ,  /*e2*/
  const int                                                       ,  /*f1*/
  const int                                                       ,  /*f2*/
  const bool                                                  /*collapsed*/
  )-> void { };
  bool ret = igl::decimate(
    VO, FO, constraint_stopping_condition,
    igl::max_faces_stopping_condition(m,orig_m,0),
    always_try,never_care,
    U, G, J, I);
  const Eigen::Array<bool,Eigen::Dynamic,1> keep = (J.array()<orig_m);
  igl::slice_mask(Eigen::MatrixXi(G),keep,1,G);
  igl::slice_mask(Eigen::VectorXi(J),keep,1,J);
  Eigen::VectorXi _1,I2;
  igl::remove_unreferenced(Eigen::MatrixXd(U),Eigen::MatrixXi(G),U,G,_1,I2);
  igl::slice(Eigen::VectorXi(I),I2,1,I);
  return ret;
}