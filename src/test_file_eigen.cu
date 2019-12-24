//****************************************************************************************
//
// Author : Aniruddha Shembekar, University of Southern California
//
//****************************************************************************************

#include "string.h"
#include <stdio.h>
#include <cmath>
#include <Eigen/Eigen>
#include "nlopt.hpp"
#include <iostream>
#include <vector>
#include "opt_obj_eigen.hpp"
#include "transformation_utilities.hpp"
#include "utilities.hpp"
#include "file_rw.hpp"
#include <chrono>

int main ()
{
    Eigen::initParallel();
    Eigen::setNbThreads(8);
    std::string file_num = "3";
    std::string file_path = "/home/aniruddha/Desktop/Composite_Layup/CUDA_experiments/data/";
    
    Eigen::MatrixXd No_particles_mat = file_rw::file_read_mat(file_path+"no_particles"+file_num+".csv");
    Eigen::MatrixXd coords = file_rw::file_read_mat(file_path+"coords"+file_num+".csv");
    Eigen::MatrixXd sheet_prop = file_rw::file_read_mat(file_path+"properties"+file_num+".csv");
    Eigen::MatrixXd node_info = file_rw::file_read_mat(file_path+"node_info"+file_num+".csv");
    Eigen::MatrixXd free_indices = file_rw::file_read_mat(file_path+"free_indices"+file_num+".csv");
    Eigen::MatrixXd X_init_mat = file_rw::file_read_mat(file_path+"X_init"+file_num+".csv");

    double OptH = 1e-11;
    double OptXtolRel = 1e-5;
    
    std::vector<double> X_init;

    for (int i=0;i<X_init_mat.cols();++i)
    {
        X_init.push_back(X_init_mat(0,i));
    }

    int No_particles = No_particles_mat(0,0);

    opt_obj::opt_obj OptObj(coords, sheet_prop, node_info, No_particles, free_indices, X_init, OptH, OptXtolRel); 
    
    std::cout << "optm start" << std::endl; 
    auto gpu_start = std::chrono::high_resolution_clock::now();
    
    bool flag = OptObj.solveOPT();
    auto gpu_end = std::chrono::high_resolution_clock::now();
    std::cout << "vector_add_gpu time: " << std::chrono::duration_cast<std::chrono::milliseconds>(gpu_end - gpu_start).count() << " milliseconds.\n";
    
    std::cout << "optm end : " << flag << std::endl; 
    
    std::vector<double> X_f_vec = OptObj.get_solx();
    Eigen::MatrixXd X_f(X_f_vec.size(),1);
    
    for (int i=0;i<X_f_vec.size();++i)
    {
        X_f(i,0) = X_f_vec[0];
    }

    double fmin = OptObj.get_solminf();
    std::cout << fmin << std::endl; 

    std::cout << "everything executed correctly..." << std::endl;

    return 0;
}