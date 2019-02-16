// #include "/usr/local/cuda-9.1/include/cuda_runtime.h"
#include "string.h"
#include <stdio.h>
#include <cmath>
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <iostream>
#include <vector>
// #include "transformation_utilities.hpp"
// #include "utilities.hpp"
// #include "file_rw.hpp"
#include <chrono>
#include "err_fun_cuda.h"


__global__ void energy_calc(double *, double *, double *, double *, double *, double, int);

double error_fun_cuda(int no_particles, Eigen::MatrixXd coords, Eigen::MatrixXd properties, Eigen::MatrixXd node_info,
    Eigen::MatrixXd free_indices, std::vector<double> x)
{
    // auto gpu_start = std::chrono::high_resolution_clock::now();
    // std::string file_path = "/home/aniruddha/Downloads/cuda_test/data/";
    // int no_particles = 160;
    // Eigen::MatrixXd coords = file_rw::file_read_mat(file_path+"coords.csv");
    // Eigen::MatrixXd properties = file_rw::file_read_mat(file_path+"properties.csv");
    // Eigen::MatrixXd node_info = file_rw::file_read_mat(file_path+"node_info.csv");
    // Eigen::MatrixXd free_indices = file_rw::file_read_mat(file_path+"free_indices.csv");
    // Eigen::MatrixXd x_mat = file_rw::file_read_mat(file_path+"X_init.csv");
    
    // std::vector<double> x;

    // for (int i=0;i<x_mat.cols();++i)
    // {
    //     x.push_back(x_mat(0,i));
    // }

    
    Eigen::MatrixXd spr_pot = Eigen::MatrixXd::Constant(no_particles,1,0);
    Eigen::MatrixXd shr_pot = Eigen::MatrixXd::Constant(no_particles,1,0);
    Eigen::MatrixXd bndg_pot = Eigen::MatrixXd::Constant(no_particles,1,0);
    int counter = -1;
    
    // for (int i = 0; i < free_indices.rows(); i++)
    // {
    //     coords(free_indices(i,0),0) = x[counter+=1];
    //     coords(free_indices(i,0),1) = x[counter+=1];
    //     coords(free_indices(i,0),2) = x[counter+=1];
    // }

    for ( int i = 0; i < free_indices.rows(); i+=3)
    {
        coords(free_indices(i,0)-1,0) = x[i];
        coords(free_indices(i,0)-1,1) = x[i+1];
        coords(free_indices(i,0)-1,2) = x[i+2];
    }

    double coords_arr[coords.rows()*coords.cols()];
    double *ptr_coords_arr;
    ptr_coords_arr = &coords_arr[0]; 
    Eigen::Map<Eigen::MatrixXd>(ptr_coords_arr,coords.rows(),coords.cols()) = coords;

    double node_info_arr[node_info.rows()*node_info.cols()];
    double *ptr_node_info_arr;
    ptr_node_info_arr = &node_info_arr[0]; 
    Eigen::Map<Eigen::MatrixXd>(ptr_node_info_arr,node_info.rows(),node_info.cols()) = node_info;    

    
    double grav_pot = coords.col(2).array().abs().sum();

    double *ptr_node_info;
    double *ptr_coords;
    double *ptr_spr_pot, *ptr_shr_pot, *ptr_bndg_pot;

    int pot_rows = spr_pot.rows();

    cudaMalloc((void**)&ptr_node_info, node_info.rows()*node_info.cols()*sizeof(double));
    cudaMalloc((void**)&ptr_coords, coords.rows()*coords.cols()*sizeof(double));
    cudaMalloc((void**)&ptr_spr_pot, pot_rows*1*sizeof(double));
    cudaMalloc((void**)&ptr_shr_pot, pot_rows*1*sizeof(double));
    cudaMalloc((void**)&ptr_bndg_pot, pot_rows*1*sizeof(double));
       
    cudaMemcpy(ptr_node_info, node_info_arr, node_info.rows()*node_info.cols()*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(ptr_coords, coords_arr, coords.rows()*coords.cols()*sizeof(double),cudaMemcpyHostToDevice);

    dim3 numBlocks(no_particles);
    dim3 threadsPerBlock(1);

    energy_calc <<<numBlocks, threadsPerBlock>>> (ptr_node_info, ptr_coords, ptr_spr_pot, ptr_shr_pot, 
        ptr_bndg_pot, properties(4,0), coords.rows());
 
    cudaMemcpy(spr_pot.data(), ptr_spr_pot, pot_rows*1*sizeof(double),cudaMemcpyDeviceToHost);
    cudaMemcpy(shr_pot.data(), ptr_shr_pot, pot_rows*1*sizeof(double),cudaMemcpyDeviceToHost);
    cudaMemcpy(bndg_pot.data(), ptr_bndg_pot, pot_rows*1*sizeof(double),cudaMemcpyDeviceToHost);
        
    cudaFree(ptr_node_info);
    cudaFree(ptr_coords);
    cudaFree(ptr_spr_pot);
    cudaFree(ptr_shr_pot);
    cudaFree(ptr_bndg_pot);

    double E = properties(0,0)*grav_pot + properties(1,0)*spr_pot.sum() + properties(2,0)*shr_pot.sum() + properties(3,0)*bndg_pot.sum();
    // std::cout << E << std::endl;
    // auto gpu_end = std::chrono::high_resolution_clock::now();
    // std::cout << "vector_add_gpu time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(gpu_end - gpu_start).count() << " nanoseconds.\n";
    return E;
}

__global__ void energy_calc(double *node_info, double *coords, 
    double *spr_pot, double *shr_pot, double *bndg_pot, double properties_val, int coords_cols)
{
    int i = blockIdx.x;
    if (i<coords_cols)
    {
        double a[3];
        double b[3];
        double c[3];
        double d[3];
        double norm_a;
        double norm_b;
        double norm_c;
        double norm_d;
        
        // Unit vectors and Spring Energy mg(dh)
        if (node_info[i]==1)
        {
            int idx = node_info[i+coords_cols*10];
            a[0] = coords[idx] - coords[i];
            a[1] = coords[idx+coords_cols] - coords[i+coords_cols];
            a[2] = coords[idx+coords_cols*2] - coords[i+coords_cols*2];    

            norm_a = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
            spr_pot[i] += (norm_a - properties_val)*(norm_a - properties_val);
            a[0] = a[0] / norm_a;
            a[1] = a[1] / norm_a;
            a[2] = a[2] / norm_a;
        }

        if (node_info[i+coords_cols]==1)
        {
            int idx = node_info[i+coords_cols*11];
            b[0] = coords[idx] - coords[i];
            b[1] = coords[idx+coords_cols] - coords[i+coords_cols];
            b[2] = coords[idx+coords_cols*2] - coords[i+coords_cols*2];
            
            norm_b = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
            spr_pot[i] += (norm_b - properties_val)*(norm_b - properties_val);
            b[0] = b[0] / norm_b;
            b[1] = b[1] / norm_b;
            b[2] = b[2] / norm_b;
        }
            
        if (node_info[i+coords_cols*2]==1)
        {
            int idx = node_info[i+coords_cols*12];
            c[0] = coords[idx] - coords[i];
            c[1] = coords[idx+coords_cols] - coords[i+coords_cols];
            c[2] = coords[idx+coords_cols*2] - coords[i+coords_cols*2];
            
            norm_c = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
            spr_pot[i] += (norm_c - properties_val)*(norm_c - properties_val);
            c[0] = c[0] / norm_c;
            c[1] = c[1] / norm_c;
            c[2] = c[2] / norm_c;
        }
        
        if (node_info[i+coords_cols*3]==1)
        {
            int idx = node_info[i+coords_cols*13];
            d[0] = coords[idx] - coords[i];
            d[1] = coords[idx+coords_cols] - coords[i+coords_cols];
            d[2] = coords[idx+coords_cols*2] - coords[i+coords_cols*2];

            norm_d = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
            spr_pot[i] += (norm_d - properties_val)*(norm_d - properties_val);
            d[0] = d[0] / norm_d;
            d[1] = d[1] / norm_d;
            d[2] = d[2] / norm_d;
        }

        // Shear Energy
        if (node_info[i+coords_cols*6]==1)
        {
            double dot_prod = a[0]*c[0] + a[1]*c[1] + a[2]*c[2];
            if (dot_prod < -1.0)
                dot_prod = -1.0;
            if (dot_prod > 1.0)
                dot_prod = 1.0;
            shr_pot[i] += tan(abs(M_PI_2 - acos(dot_prod)));
        }

        if (node_info[i+coords_cols*7]==1)
        {
            double dot_prod = b[0]*c[0] + b[1]*c[1] + b[2]*c[2];
            if (dot_prod < -1.0)
                dot_prod = -1.0;
            if (dot_prod > 1.0)
                dot_prod = 1.0;
            shr_pot[i] += tan(abs(M_PI_2 - acos(dot_prod)));
        }
        
        if (node_info[i+coords_cols*8]==1)
        {
            double dot_prod = b[0]*d[0] + b[1]*d[1] + b[2]*d[2];
            if (dot_prod < -1.0)
                dot_prod = -1.0;
            if (dot_prod > 1.0)
                dot_prod = 1.0;
            shr_pot[i] += tan(abs(M_PI_2 - acos(dot_prod)));
        }
        
        if (node_info[i+coords_cols*9]==1)
        {
            double dot_prod = a[0]*d[0] + a[1]*d[1] + a[2]*d[2];
            if (dot_prod < -1.0)
                dot_prod = -1.0;
            if (dot_prod > 1.0)
                dot_prod = 1.0;
            shr_pot[i] += tan(abs(M_PI_2 - acos(dot_prod)));
        }

        // Bending Energy
        if (node_info[i+coords_cols*4]==1)
        {
            double dot_prod = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
            if (dot_prod < -1.0)
                dot_prod = -1.0;
            if (dot_prod > 1.0)
                dot_prod = 1.0;
            bndg_pot[i] += tan((M_PI - acos(dot_prod))/2);
        }
        
        if (node_info[i+coords_cols*5]==1)
        {
            double dot_prod = d[0]*c[0] + d[1]*c[1] + d[2]*c[2];
            if (dot_prod < -1.0)
                dot_prod = -1.0;
            if (dot_prod > 1.0)
                dot_prod = 1.0;
            bndg_pot[i] += tan((M_PI - acos(dot_prod))/2);
        }
    }
}
