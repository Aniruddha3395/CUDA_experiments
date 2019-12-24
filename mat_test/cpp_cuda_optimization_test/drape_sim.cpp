// Author   : Aniruddha Shembekar, Research Engineer, University of Southern California

// #include "compute_energy.hpp"
#include "mex.h"
#include "matrix.h"
#include "string.h"
#include <stdio.h>
#include <cmath>
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <nlopt.hpp>
#include <iostream>
#include <vector>
#include <iterator>
#include <chrono>
#include "opt_obj.h"

void mexFunction (int _no_op_args, mxArray *mex_op[], int _no_ip_args, const mxArray *mex_in[] )
{
    /////////////////// DECLARATIONS ///////////////////////////////////
    // I/O variables

    /////////////////// INPUT ///////////////////////////////////
    double* ptr_no_particles;
    ptr_no_particles = mxGetDoubles(mex_in[3]);
    int No_particles = (int) *ptr_no_particles;
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> coords (mxGetPr(mex_in[0]), No_particles, 3);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> sheet_prop (mxGetPr(mex_in[1]), 5, 1);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> node_info (mxGetPr(mex_in[2]), No_particles, 14);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> free_indices (mxGetPr(mex_in[4]), mxGetNumberOfElements(mex_in[4]), 1);
    double* X_init;
    X_init = mxGetDoubles(mex_in[5]);
    std::vector<double> X_init_vec(X_init, X_init+mxGetNumberOfElements(mex_in[5]));
    
    /////////////////// CREATING OPTIMIZATION OBJECT //////////////////

    double OptH = 1e-11;
    double OptXtolRel = 1e-5;
    
    double *ptr_node_info_arr;
    double node_info_arr[node_info.rows()*node_info.cols()];
    ptr_node_info_arr = &node_info_arr[0]; 
    Eigen::Map<Eigen::MatrixXd>(ptr_node_info_arr,node_info.rows(),node_info.cols()) = node_info;   

    int node_info_size = node_info.rows()*node_info.cols();

    opt_obj::opt_obj OptObj(coords, sheet_prop, No_particles, free_indices, X_init_vec, OptH, OptXtolRel, node_info_arr, node_info_size);
    
//     auto gpu_start = std::chrono::high_resolution_clock::now();
    bool flag = OptObj.solveOPT();
//     auto gpu_end = std::chrono::high_resolution_clock::now();
//     std::cout << "vector_add_gpu time: " << std::chrono::duration_cast<std::chrono::milliseconds>(gpu_end - gpu_start).count() << " milliseconds.\n";
        
    std::vector<double> X_f_vec = OptObj.get_solx();
    Eigen::MatrixXd X_f(X_f_vec.size(),1);
    
    for (int i=0;i<X_f_vec.size();++i)
    {
        X_f(i,0) = X_f_vec[i];
    }

    double fmin = OptObj.get_solminf();

    /////////////////// OUTPUT ///////////////////////////////////
    mex_op[0] = mxCreateDoubleMatrix(X_f.rows(),X_f.cols(), mxREAL);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M0 (mxGetPr(mex_op[0]),X_f.rows(),X_f.cols());
    M0 = X_f;
    mex_op[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M1 (mxGetPr(mex_op[1]),1,1);
    M1 << fmin;  
    return;
}
