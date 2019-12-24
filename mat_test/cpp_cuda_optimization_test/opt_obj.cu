// Author   : Aniruddha Shembekar, Research Engineer, University of Southern California

#include <iostream>
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <vector>
#include <string>
#include <cmath>
#include <nlopt.hpp>
#include "stdlib.h"
#include "opt_obj.h"
#include <chrono>

__global__ void energy_calc(double *node_info, double *coords, 
    double *spr_pot, double *shr_pot, double *bndg_pot, double properties_val, int coords_cols)
{
    int i = threadIdx.x + 1000*blockIdx.x;
    
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
        spr_pot[i] = 0;
        shr_pot[i] = 0;
        bndg_pot[i] = 0;
        double angle;

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
//             bndg_pot[i] += tan((M_PI - acos(dot_prod))/2);
            angle = acos(dot_prod);
            bndg_pot[i] += 1.124e10 * exp(-15.34 * angle) + 0.06109 * exp(-1.1386 * angle);
        }
        
        if (node_info[i+coords_cols*5]==1)
        {
            double dot_prod = d[0]*c[0] + d[1]*c[1] + d[2]*c[2];
            if (dot_prod < -1.0)
                dot_prod = -1.0;
            if (dot_prod > 1.0)
                dot_prod = 1.0;
//             bndg_pot[i] += tan((M_PI - acos(dot_prod))/2);
            angle = acos(dot_prod);
            bndg_pot[i] += 1.124e10 * exp(-15.34 * angle) + 0.06109 * exp(-1.1386 * angle);
        }
    }

}


namespace opt_obj
{
    // defualt nlopt fuction
    double customminfunc(const std::vector<double>& x, std::vector<double>& grad, void* data) {
        // Because we wanted a Class
        // without static members, but NLOpt library does not support
        // passing methods of Classes, we use these auxilary functions.
        opt_obj *c = (opt_obj *) data;
        return c->ObjFun(x,grad);
    }

    // class constructor
    opt_obj::opt_obj( Eigen::MatrixXd _coords, Eigen::MatrixXd _prop, int _NoParticles, Eigen::MatrixXd _freeIndices, 
        std::vector<double> _X_init, double OptH, double OptXtolRel, double node_info_arr[], int _node_info_size) 
    {

        //choose optimizer
        // alg_type = nlopt::LN_NEWUOA;
        // optalg = nlopt::LN_NEWUOA_BOUND;
        // alg_type = nlopt::LN_BOBYQA;
        // alg_type = nlopt::LN_COBYLA;
        // alg_type = nlopt::LD_SLSQP;
        alg_type = nlopt::LD_LBFGS;
        // optalg = nlopt::GN_ISRES;
        
        
        coords = _coords;
        properties = _prop;
        node_info_size = _node_info_size;
        no_particles = _NoParticles;
        free_indices = _freeIndices;
        X_init = _X_init;
        optH = OptH;
        optXtolRel = OptXtolRel;
        solminf = 0;
        ptr_node_info_arr = &node_info_arr[0];

        // Optimization Parameters
        OptVarDim = _X_init.size();
        opt = nlopt::opt(alg_type, OptVarDim);
        // OptVarlb.resize(OptVarDim);
        // OptVarub.resize(OptVarDim);
        
        // for (int i=0;i<OptVarlb.size();++i)
        // {
        //     OptVarlb[i] = -100000;
        //     OptVarub[i] = 100000;
        // }

        opt.set_xtol_rel(optXtolRel);
        opt.set_min_objective(customminfunc, this);

        // opt.set_lower_bounds(OptVarlb);
        // opt.set_upper_bounds(OptVarub);
        opt.set_maxeval(50);
        
        // optimization params
        // OptVarDim = x_start.size();

        
        // opt.set_xtol_rel(optXtolRel);
        optx.resize(OptVarDim);
        // opt.set_min_objective(customminfunc, this);
    }

    // class destructor
    opt_obj::~opt_obj()
    {
    }

    
    double opt_obj::ErrFun_cuda(const std::vector<double> &x)
    {
        int counter = -1;
        for ( int i = 0; i < free_indices.rows(); i++)
        {
            coords(free_indices(i,0),0) = x[counter+=1];
            coords(free_indices(i,0),1) = x[counter+=1];
            coords(free_indices(i,0),2) = x[counter+=1];
        }
        
        
        cudaMemcpy(ptr_coords, coords.data(), coords.rows()*coords.cols()*sizeof(double),cudaMemcpyHostToDevice);
        
        energy_calc <<< block_count+1, no_particles%1000 >>> (ptr_node_info, ptr_coords, ptr_spr_pot, ptr_shr_pot, ptr_bndg_pot, properties(4,0), coords.rows());
        
        
        // auto gpu_start = std::chrono::high_resolution_clock::now();
        cudaMemcpy(shr_pot.data(), ptr_shr_pot, no_particles*sizeof(double),cudaMemcpyDeviceToHost);
        cudaMemcpy(spr_pot.data(), ptr_spr_pot, no_particles*sizeof(double),cudaMemcpyDeviceToHost);
        cudaMemcpy(bndg_pot.data(), ptr_bndg_pot, no_particles*sizeof(double),cudaMemcpyDeviceToHost);
        // auto gpu_end = std::chrono::high_resolution_clock::now();
        // std::cout << " time 2 : " << std::chrono::duration_cast<std::chrono::microseconds>(gpu_end - gpu_start).count() << " microseconds.\n";
        
        std::cout << properties(0,0)*coords.col(2).array().abs().sum() + properties(1,0)*spr_pot.sum() + properties(2,0)*shr_pot.sum() + properties(3,0)*bndg_pot.sum() << std::endl;
        return properties(0,0)*coords.col(2).array().abs().sum() + properties(1,0)*spr_pot.sum() + properties(2,0)*shr_pot.sum() + properties(3,0)*bndg_pot.sum();
        
    }

    bool opt_obj::solveOPT()
    {
        spr_pot = Eigen::MatrixXd::Constant(no_particles,1,0);
        shr_pot = Eigen::MatrixXd::Constant(no_particles,1,0);
        bndg_pot = Eigen::MatrixXd::Constant(no_particles,1,0);
        
        cudaMalloc((void**)&ptr_spr_pot, no_particles*1*sizeof(double));
        cudaMalloc((void**)&ptr_shr_pot, no_particles*1*sizeof(double));
        cudaMalloc((void**)&ptr_bndg_pot, no_particles*1*sizeof(double));
        cudaMalloc((void**)&ptr_coords, coords.rows()*coords.cols()*sizeof(double));
        
        cudaMalloc((void**)&ptr_node_info, node_info_size*sizeof(double));
        cudaMemcpy(ptr_node_info, ptr_node_info_arr, node_info_size*sizeof(double),cudaMemcpyHostToDevice);
        
        block_count = no_particles/1000;
        solx = X_init;
        bool successFlag = false;
        try
        {    
            nlopt::result result = opt.optimize(solx, solminf);
            successFlag = true;
        }
        catch(std::exception &e) {
            std::cout << "nlopt failed: " << e.what() << std::endl;
        }

        cudaFree(ptr_node_info);
        cudaFree(ptr_spr_pot);
        cudaFree(ptr_shr_pot);
        cudaFree(ptr_bndg_pot);
        cudaFree(ptr_coords);
        
        return successFlag;
    };

    std::vector<double> opt_obj::get_solx()
    {
        return solx;
    };

    double opt_obj::get_solminf()
    {
        return solminf;
    };

    // gradient computation:
    // Forward Difference Method
    double opt_obj::ObjFun(const std::vector<double> &x, std::vector<double> &grad)
    {
        double err = ErrFun_cuda(x);
        if (!grad.empty()) {
            std::vector<double> xph = x;
            for (uint i=0; i < x.size(); ++i)
            {
                xph[i] += optH;
                grad[i] = (ErrFun_cuda(xph)-err)/optH;
                xph[i] -= optH;
            }
        }    
        return err;
    };
}

