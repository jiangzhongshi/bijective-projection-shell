#include "util.h"
#include <igl/boundary_loop.h>
#include <igl/list_to_matrix.h>
#include <map>

// Code from Hanxiao

void split_along_edges(
    int N,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    const Eigen::VectorXi& B,
    int si,
    int ti,
    Eigen::VectorXi& nB
){
    int num = V.rows();
    int numf = F.rows();
    long n_f = F.rows();
    long n_v = V.rows();
    long dim = V.cols();
    std::vector<int> nBvec;
    int sz;
    // number of edges available
    sz = si < ti ? ti-si: ti+B.rows()-si;
    std::vector<int> numbers(sz);
    double total = 0;
    int addup = 0;
    double e = double(N)/sz;
    for(int i=0;i<sz;i++){
        numbers[i] = int(total + e)-addup; 
        total += e;
        addup += numbers[i];
    }
    if(addup!=N)
        numbers.back() = N-(addup-numbers.back());
    V.conservativeResize(V.rows()+2*N,dim);
    F.conservativeResize(F.rows()+2*(N+sz),3);
    // remember the points coord for pts already inserted
    std::map<std::pair<int,int>,Eigen::MatrixXd> points_added;
    std::map<std::pair<int,int>,std::vector<int>> points_added_i;
    for(int i=0;i<sz;i++){
        int a = B((si+i)%B.rows());
        int b = B((si+i+1)%B.rows());
        int t1,t2,t3;
        int n = numbers[i];
        if(n==0)continue;
        for(int j=0;j<n_f;j++){
            int index = -1;
            for(int k=0;k<3;k++){
                if((F(j,k)==a && F(j,(k+1)%3)==b) ||
                   (F(j,k)==b && F(j,(k+1)%3)==a)){
                       index = k;
                       break;
                }
            }
            if(index == -1) continue;
            t1 = F(j,index);
            t2 = F(j,(index+1)%3);
            t3 = F.row(j).sum()-t1-t2;
            int start = n_v;
            std::pair<int,int> c(t2,t1);
            std::vector<int> vl;
            if(points_added.find(c)==points_added.end()){
                c = std::make_pair(t1,t2);
                points_added[c]=Eigen::MatrixXd(n,dim);
                for(int k=0;k<n;k++){
                    V.row(n_v)=(V.row(t2)-V.row(t1)).array()*(k+1)/(n+1)+V.row(t1).array();
                    points_added[c].row(k)=V.row(n_v);
                    points_added_i[c].push_back(n_v);
                    n_v++;
                }
                vl.push_back(t1);
                vl.insert(vl.end(),points_added_i[c].begin(),points_added_i[c].end());
                vl.push_back(t2);
            }else{
                vl.push_back(t1);
                std::vector<int> tmp = points_added_i[c];
                std::reverse(tmp.begin(),tmp.end());
                vl.insert(vl.end(),tmp.begin(),tmp.end());
                vl.push_back(t2);                 
            }
            for(int k=0;k<n+1;k++){
                F.row(n_f++)<<t3,vl[k],vl[k+1];
            }
            F.block(j,0,F.rows()-j-1,3) = F.bottomRows(F.rows()-j-1);
            F.conservativeResize(F.rows()-1,3);
            n_f-=1;
            j-=1;
        }
    }
    for(int i=0;i<sz;i++){
        nBvec.push_back(B((si+i)%B.rows()));
        std::pair<int,int> c1(B((si+i)%B.rows()),B((si+i+1)%B.rows()));
        std::pair<int,int> c2(B((si+i+1)%B.rows()),B((si+i)%B.rows()));
        if(points_added_i.find(c1)!=points_added_i.end())
            nBvec.insert(nBvec.end(),points_added_i[c1].begin(),points_added_i[c1].end());
        else{
            auto tmp = points_added_i[c2];
            std::reverse(tmp.begin(),tmp.end());
            nBvec.insert(nBvec.end(),tmp.begin(),tmp.end());
        }
    }
    nBvec.push_back(B(B.rows()-1));
    igl::list_to_matrix(nBvec,nB);
    V.conservativeResize(n_v,dim);
    F.conservativeResize(n_f,3);
}



// split edges on the boundary
void split_boundary(
    const std::vector<int>& n,
    const std::vector<std::pair<int,int>>& P,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
){
    std::vector<std::vector<int>> Bv;
    igl::boundary_loop(F,Bv);
    for(int p=0;p<P.size();p++){
        int a = P[p].first;
        int b = P[p].second;
        int pos_a=-1,pos_b=-1;
        int bdi=-1;
        for(int i=0;i<Bv.size();i++){
            for(int j=0;j<Bv[i].size();j++){
                if(Bv[i][j]==a){
                    pos_a = j;
                }
                if(Bv[i][j]==b){
                    pos_b = j;
                }
            }
            if(pos_a!=-1&&pos_b!=-1){
                bdi = i;
                break;
            }
        }
        assert(bdi!=-1 && "did not match boundary");
        Eigen::VectorXi B,_nB;
        igl::list_to_matrix(Bv[bdi],B);
        split_along_edges(n[p],V,F,B,pos_a,pos_b,_nB);
        Bv.clear();
        igl::boundary_loop(F,Bv);
    }
}

#include <igl/flip_avoiding_line_search.h>
namespace igl::flip_avoiding{
    double compute_max_step_from_singularities(const Eigen::MatrixXd& uv,
                                                          const Eigen::MatrixXi& F,
                                                          Eigen::MatrixXd& d);
}
std::tuple<double, double, bool> flip_avoiding_line_search(
  const Eigen::MatrixXi F,
  Eigen::MatrixXd& cur_v,
  Eigen::MatrixXd& dst_v,
  int max_iter,
  std::function<double(Eigen::MatrixXd&)> energy,
  double cur_energy)
{
  // This function is adapted from libigl

  Eigen::MatrixXd d = dst_v - cur_v;
    
  double min_step_to_singularity = igl::flip_avoiding::compute_max_step_from_singularities(cur_v,F,d);
  double step_size = std::min(1., min_step_to_singularity*0.8);

  double old_energy;
  if (cur_energy > 0)
  {
    old_energy = cur_energy;
  }
  else
  {
    old_energy = energy(cur_v); // no energy was given -> need to compute the current energy
  }
  double new_energy = old_energy;
  int cur_iter = 0; int MAX_STEP_SIZE_ITER = max_iter;

  while (new_energy >= old_energy && cur_iter < MAX_STEP_SIZE_ITER)
  {
    Eigen::MatrixXd new_x = cur_v + step_size * d;

    double cur_e = energy(new_x);
    if ( cur_e >= old_energy)
    {
      step_size /= 2;
    }
    else
    {
      cur_v = new_x;
      new_energy = cur_e;
    }
    cur_iter++;
  }
  return std::make_tuple(new_energy, step_size, cur_iter!=max_iter);
}