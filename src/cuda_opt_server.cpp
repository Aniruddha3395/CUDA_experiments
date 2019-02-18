#include "ros/ros.h"
#include "cuda_opt/OptCuda.h"
#include "string.h"
#include <stdio.h>
#include <cmath>
#include <Eigen/Eigen>
#include "nlopt.hpp"
#include <iostream>
#include <vector>
#include "opt_obj.h"
#include "utilities.hpp"
#include "file_rw.hpp"
#include <chrono>

bool cuda_opt_fun(cuda_opt::OptCuda::Request &req, cuda_opt::OptCuda::Response &res)
{
    std::string file_num = req.file_num;
    std::string file_path = "/home/aniruddha/Desktop/Composite_Layup/CUDA_experiments/mat_test/cpp_cuda_optimization_test/data_files/";
    
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
    int node_info_size= node_info.rows()*node_info.cols();

    double *ptr_node_info_arr;
    double node_info_arr[node_info.rows()*node_info.cols()];
    ptr_node_info_arr = &node_info_arr[0]; 
    Eigen::Map<Eigen::MatrixXd>(ptr_node_info_arr,node_info.rows(),node_info.cols()) = node_info;   

    opt_obj::opt_obj OptObj(coords, sheet_prop, No_particles, free_indices, X_init, OptH, OptXtolRel, node_info_arr, node_info_size); 
    
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

    res.error = OptObj.get_solminf();

    // res.error = 111;
    std::cout << res.error << std::endl; 

    std::cout << "everything executed correctly..." << std::endl;

    return true;
}

int main(int argc, char **argv)
{

    ros::init(argc, argv, "drape_sim_server");
    ros::NodeHandle n;
    ros::ServiceServer service = n.advertiseService("opt_w_cuda", cuda_opt_fun);
    ROS_INFO("Ready to do some optimization");
    ros::spin();


    return 0;
}
