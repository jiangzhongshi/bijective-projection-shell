#include <igl/matrix_to_list.h>
#include <igl/serialize.h>
#include "matchmaker.h"
#include "edge_split.h"
#include "decompose_polygon.h"
#include "mst.h"
#include "cut_mesh/HalfEdgeIterator.h"

#include "target_polygon.h"
#include "progressive_embedding.h"
#include "loader.h"
#include "argh.h"
#include "cut_mesh/edge_flaps.h"
#include "local_smooth/local_smooth.h"
#include <igl/remove_unreferenced.h>
#include <igl/copyleft/cgal/orient2D.h>
#include "decompose_polygon.h"
#include "mst.h"
#include "cut_mesh/HalfEdgeIterator.h"
#include "edge_split.h"

#include <igl/vertex_triangle_adjacency.h>
#include <igl/boundary_loop.h>
#include <igl/matrix_to_list.h>
#include <igl/doublearea.h>

#include <unordered_map>
#include <igl/writeDMAT.h>
#include <igl/writeOFF.h>
#include <igl/writeOBJ.h>
#include <algorithm>
#include <unordered_set>
#include "progressive_embedding.h"
#include "validity_check.h"
#include <igl/slim.h>

// based on the polygon list build a graph
extern void build_graph(
    int n, // size of boundary
    const Eigen::MatrixXd& V,
    const std::vector<std::vector<int>>& L,
    Eigen::SparseMatrix<int>& graph);

extern void remove_ears(
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F
    );

extern void neighbor_sector_bar(
        const std::vector<std::vector<int>>& L, // polygon list
        const Eigen::MatrixXd& V2,
        const Eigen::MatrixXi& F2,
        const std::set<std::pair<int,int>>& sector_bar,
        int s,
        int t,
        std::pair<int,int>& b00,
        std::pair<int,int>& b01
    );

extern bool match_sector_bar(
        std::map<std::pair<int,int>,std::vector<int>>& splits,
        const Eigen::VectorXi& T,
        std::pair<int,int>& b,
        std::pair<int,int>& c,
        int center
    );

extern void collect_blocked_face(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        std::set<int>& no_enter_f,
        const std::pair<int,int>& c00,
        const std::pair<int,int>& c01
    );

extern void mark_impassible(
        const std::vector<std::vector<int>>& L,
        const Eigen::MatrixXd& V2,
        const Eigen::MatrixXi& F2,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::VectorXi& T, // match 2d to 3d
        const std::set<std::pair<int,int>>& sector_bar,
        const std::pair<int,int>& t, // trace [a,b]
        std::map<std::pair<int,int>,std::vector<int>>& splits, // path traced in 3d
        std::set<int>& no_enter_f
    );

extern void prapare_TT(
        const Eigen::MatrixXi& F,
        const std::set<int>& no_enter_f,
        const std::vector<std::pair<int,int>>& M,
        const std::vector<std::vector<int>>& VF,
        Eigen::MatrixXi& TT,
        Eigen::MatrixXi& TTi
    );

void match_to_square_with_constraints(
  Eigen::MatrixXd& V, // 3d mesh, along with result
  Eigen::MatrixXi& F, // 3d mesh, along with result
  const Eigen::MatrixXi& F2,  // 2d coarse mesh (has to contain boundary)
  const std::vector<std::vector<int>>& L, // polygon
  const Eigen::VectorXi & cons_i, //
  const std::vector<std::pair<int,int>>& to_trace, // interior edges, to be traced
  Eigen::VectorXi& bound, // input: boundary indices, return: cut
  std::vector<std::vector<int>> &total_splits
){

  Eigen::VectorXi T;
  igl::cat(1, bound, cons_i, T);

  int bod_num = bound.size();
    // update corresponding map

    // [ by the merit of shor algorithm, we know the ]
    // [ range of boundary vertices are 0: bod_num-1 ]

    std::set<int> no_enter; // set of vertices that should avoid
    std::vector<std::pair<int,int>> mask_e;
    std::map<std::pair<int,int>,std::vector<int>> splits;
    Eigen::VectorXi bi;
    igl::boundary_loop(F,bi);
    for(int i=0;i<bi.rows();i++)
        no_enter.insert(bi(i));

    // [ traces already in the mesh together with
    //   boundary edges forms a set of sectors ]
    std::set<std::pair<int,int>> sector_bar;
    for(auto e: to_trace){
        for(int i : {e.first,e.second}){
            if(i<bod_num){
                int i_prev = (i-1+bod_num)%bod_num;
                int i_next = (i+1)%bod_num;
                sector_bar.insert(std::make_pair(i,i_prev));
                sector_bar.insert(std::make_pair(i_next,i));
            }
        }
    }

    std::vector<std::pair<int,int>> to_trace_ordered(to_trace.begin(), to_trace.end());

    std::set<std::pair<int,int>> on_path;
    std::set<int> on_path_v;
    for(int i=0;i<bound.rows();i++){
        int a = std::min(bound(i),bound((i+1)%bound.rows()));
        int b = std::max(bound(i),bound((i+1)%bound.rows()));
        on_path.insert(std::make_pair(a,b));
        on_path_v.insert(bound(i));
    }
    for (int i=0; i<cons_i.rows(); i++) {
        on_path_v.insert(cons_i[i]); // Bugfix by Zhongshi at May 30. split when the cons is on neck
    }
    Eigen::MatrixXi TT,TTi;
    igl::triangle_triangle_adjacency(F,TT,TTi);
    std::vector<std::vector<int>> VF,VFi;
    igl::vertex_triangle_adjacency(V,F,VF,VFi);
    int xc = 0;
    for(int i=0;i<cons_i.rows();i++){
        no_enter.insert(cons_i(i));
    }
    for(auto e: to_trace_ordered){
        int l = e.first;
        int r = e.second;
        std::cout<<"traced 2d "<<l<<" to "<<r<<std::endl;
        std::cout<<"traced 3d "<<T(l)<<" to "<<T(r)<<std::endl;
        // start and end point should not be avoided
        no_enter.erase(T(l));
        no_enter.erase(T(r));

        std::set<int> no_enter_f;
        Eigen::MatrixXd V2;
        mark_impassible(L,V2,F2,V,F,T,sector_bar,e,splits,no_enter_f);
        auto Vn = V;
        auto Fn = F;
        std::vector<int> E; // traced path
        bool trace_succeed = path_tracing(V,F,std::make_pair(T(l),T(r)),no_enter,mask_e,no_enter_f,TT,Vn,Fn,E);
        no_enter.insert(T(l));
        no_enter.insert(T(r));

        if (!trace_succeed) continue;
        std::reverse(E.begin(),E.end());
        splits[e] = E;
        assert(E.size()>=2);
        assert(T(e.first) == E[0] && T(e.second) == E.back());


        // update impassible edges info
        for(int i=0;i<E.size();i++){
            no_enter.insert(E[i]);
            if(i == E.size()-1) break;
            int a = E[i];
            int b = E[i+1];
            mask_e.push_back(std::make_pair(a,b));
            on_path.insert(std::make_pair(std::min(a,b),std::max(a,b)));
            on_path_v.insert(a);
            on_path_v.insert(b);
        }
        sector_bar.insert(e);
        //igl::triangle_triangle_adjacency(F,TT,TTi);
        igl::vertex_triangle_adjacency(V,F,VF,VFi);
        no_enter_f.clear();
        prapare_TT(F,no_enter_f,mask_e,VF,TT,TTi);
        // split edge whose both end points are on boundary
        std::vector<std::vector<double>> V_vec;
        std::vector<std::vector<int>> F_vec;
        std::vector<std::vector<int>> TT_vec;
        std::vector<std::vector<int>> TTi_vec;
        igl::matrix_to_list(V,V_vec);
        igl::matrix_to_list(F,F_vec);
        igl::matrix_to_list(TT,TT_vec);
        igl::matrix_to_list(TTi,TTi_vec);
        for(int i=0;i<F_vec.size();i++){
            for(int k=0;k<3;k++){
                int a1 = std::min(F_vec[i][k],F_vec[i][(k+1)%3]);
                int a2 = std::max(F_vec[i][k],F_vec[i][(k+1)%3]);
                auto t = std::make_pair(a1,a2);
                if(on_path.find(t) != on_path.end()) continue;
                if(on_path_v.find(F_vec[i][k])!=on_path_v.end() &&
                   on_path_v.find(F_vec[i][(k+1)%3])!=on_path_v.end()){
                    bool status = edge_split(V_vec,F_vec,TT_vec,TTi_vec,i,k);
                    assert(status!=false && "split boundary");
                }
            }
        }
        igl::list_to_matrix(TT_vec,TT);
        igl::list_to_matrix(TTi_vec,TTi);
        igl::list_to_matrix(F_vec,F);
        igl::list_to_matrix(V_vec,V);
        std::cout<<"iteration "<<xc++<<"/"<<to_trace_ordered.size()<<std::endl;
        std::cout<<"current #F "<<F.rows()<<std::endl;
    }

    // sample on to_trace
    Eigen::MatrixXd known;
    std::vector<std::vector<double>> known_vec;
    Eigen::VectorXi ki;
    std::vector<int> ki_vec;
    for(auto e: to_trace_ordered){
        int a = e.first;
        int b = e.second;
        // on 2d [a     -     b]
        // on 3d [T(a) --- T(b)]
        std::vector<int> path = splits[e];
        total_splits.push_back(splits[e]);
            
        for(int k=0;k<path.size()-2;k++){
            ki_vec.push_back(path[k+1]);
        }
    }
    igl::list_to_matrix(ki_vec,ki);
    T.conservativeResize(T.rows()+ki.rows());
    T.bottomRows(ki.rows())<<ki;
    bound = T;
}
