#include "compute_energy.hpp"
#include "mex.h"
#include "matrix.h"
#include "string.h"
#include <stdio.h>
#include <cmath>
#include </usr/local/include/eigen3/Eigen/Eigen>
#include </usr/local/include/nlopt.h>
#include <iostream>
#include <vector>
#include "nlopt.hpp"

// default nlopt fuction
double customminfunc( const std::vector<double>& x, std::vector<double>& grad, void* data ) {
    // Because we wanted a Class
    // without static members, but NLOpt library does not support
    // passing methods of Classes, we use these auxilary functions.
    opt_obj *c = (opt_obj *) data;
    return c->ObjFun(x,grad);
};
    
    
class opt_obj
{
    int OptVarDim;
    std::vector<double> optx;
    std::vector<double> solx;
    double solminf;
    nlopt::algorithm alg_type;
    nlopt::opt opt;
    Eigen::MatrixXd coords;
    Eigen::VectorXd properties;
    Eigen::MatrixXd node_info;
    int no_particles;
    Eigen::VectorXd free_indices;
    Eigen::VectorXd X_init;
    
    // Forward Difference Method
    double ObjFun(const std::vector<double> &x, std::vector<double> &grad)
    {
        double err = ErrFun(x);
        if (!grad.empty()) {
            std::vector<double> xph = x;
            for (uint i=0; i < x.size(); ++i)
            {
                xph[i] += 1e-9;
                grad[i] = (ErrFun(xph)-err)/1e-9;
                xph[i] -= 1e-9;
            }
        }
        return err;
    };
    
    opt_obj( Eigen::MatrixXd _coords, Eigen::VectorXd _prop, Eigen::MatrixXd _nodeInfo, int _NoParticles, Eigen::VectorXd _freeIndices, Eigen::VectorXd _X_init )
    {
        //choose optimizer
        // alg_type = nlopt::LN_NEWUOA;
        // alg_type = nlopt::LN_NEWUOA_BOUND;
        // alg_type = nlopt::LN_BOBYQA;
        // alg_type = nlopt::LN_COBYLA;
        alg_type = nlopt::LD_LBFGS;
        // alg_type = nlopt::GN_ISRES;
        
        
        coords = _coords;
        properties = _prop;
        node_info = _nodeInfo;
        no_particles = _NoParticles;
        free_indices = _freeIndices;
        X_init = _X_init;
        
        OptVarDim = _X_init.size();
        opt = nlopt::opt(alg_type, OptVarDim);
        
        // Optimization Parameters
        opt.set_xtol_rel( 1e-6 );
        opt.set_min_objective(customminfunc, this);
        
    };
    
    // class destructor
    ~opt_obj()
    {
    };
    
    
    double ErrFun(const std::vector<double> &x)
    {
        double grav_pot = 0;
        double spr_pot = 0;
        double shr_pot = 0;
        double bndg_pot = 0;
        Eigen::MatrixXd a(1,3);
        Eigen::MatrixXd b(1,3);
        Eigen::MatrixXd c(1,3);
        Eigen::MatrixXd d(1,3);
        double norm_a;
        double norm_b;
        double norm_c;
        double norm_d;
        double dot_prod;
        
        for ( int i = 0; i < free_indices.size(); i+3 )
        {
            coords(free_indices(i),0) = x[i];
            coords(free_indices(i),1) = x[i+1];
            coords(free_indices(i),2) = x[i+2];
        }
        
        
        /////////////////// DEFINITION ///////////////////////////////////
        for ( int i = 0; i < no_particles; i++ )
        {
            grav_pot += std::abs( coords(i,2) );
            
            // Unit vectors and Spring Energy mg(dh)
            if (node_info(i,0)==1)
            {
                a << coords.row(node_info(i,10)).array() - coords.row(i).array();
                norm_a = a.norm();
                spr_pot += pow( norm_a - properties(4,0), 2);
                a = a / norm_a;
            }
            
            if (node_info(i,1)==1)
            {
                b << coords.row(node_info(i,11)).array() - coords.row(i).array();
                norm_b = b.norm();
                spr_pot += pow( norm_b - properties(4,0), 2);
                b = b / norm_b;
            }
            
            if (node_info(i,2)==1)
            {
                c << coords.row(node_info(i,12)).array() - coords.row(i).array();
                norm_c = c.norm();
                spr_pot += pow( norm_c - properties(4,0), 2);
                c = c / norm_c;
            }
            
            if (node_info(i,3)==1)
            {
                d << coords.row(node_info(i,13)).array() - coords.row(i).array();
                norm_d = d.norm();
                spr_pot += pow( norm_d - properties(4,0), 2);
                d = d / norm_d;
            }
            
            // Shear Energy
            if (node_info(i,6)==1)
            {
                dot_prod = a(0,0)*c(0,0) + a(0,1)*c(0,1) + a(0,2)*c(0,2);
                if (dot_prod < -1.0)
                    dot_prod = -1.0;
                if (dot_prod > 1.0)
                    dot_prod = 1.0;
                shr_pot += tan( std::abs(M_PI_2 - acos( dot_prod )) );
            }
            
            if (node_info(i,7)==1)
            {
                dot_prod = b(0,0)*c(0,0) + b(0,1)*c(0,1) + b(0,2)*c(0,2);
                if (dot_prod < -1.0)
                    dot_prod = -1.0;
                if (dot_prod > 1.0)
                    dot_prod = 1.0;
                shr_pot += tan( std::abs(M_PI_2 - acos( dot_prod )) );
            }
            
            if (node_info(i,8)==1)
            {
                dot_prod = b(0,0)*d(0,0) + b(0,1)*d(0,1) + b(0,2)*d(0,2);
                if (dot_prod < -1.0)
                    dot_prod = -1.0;
                if (dot_prod > 1.0)
                    dot_prod = 1.0;
                shr_pot += tan( std::abs(M_PI_2 - acos( dot_prod )) );
            }
            
            if (node_info(i,9)==1)
            {
                dot_prod = a(0,0)*d(0,0) + a(0,1)*d(0,1) + a(0,2)*d(0,2);
                if (dot_prod < -1.0)
                    dot_prod = -1.0;
                if (dot_prod > 1.0)
                    dot_prod = 1.0;
                shr_pot += tan( std::abs(M_PI_2 - acos( dot_prod )) );
            }
            
            // Bending Energy
            if (node_info(i,4)==1)
            {
                dot_prod = a(0,0)*b(0,0) + a(0,1)*b(0,1) + a(0,2)*b(0,2);
                if (dot_prod < -1.0)
                    dot_prod = -1.0;
                if (dot_prod > 1.0)
                    dot_prod = 1.0;
                bndg_pot += tan( (M_PI - acos( dot_prod ))/2 );
            }
            
            if (node_info(i,5)==1)
            {
                dot_prod = d(0,0)*c(0,0) + d(0,1)*c(0,1) + d(0,2)*c(0,2);
                if (dot_prod < -1.0)
                    dot_prod = -1.0;
                if (dot_prod > 1.0)
                    dot_prod = 1.0;
                bndg_pot += tan( (M_PI - acos( dot_prod ) )/2 );
            }
        }
        
        return properties(0,0) * grav_pot + properties(1,0) * spr_pot + properties(2,0) * shr_pot + properties(3,0) * bndg_pot;
    };
    
    bool solveOPT()
    {
        solx = X_init;
        bool successFlag = false;
        try{
            nlopt::result result = opt.optimize(solx, solminf);
            successFlag = true;
        }
        catch(std::exception &e) {
            std::cout << "nlopt failed: " << e.what() << std::endl;
        }
        return successFlag;
    };
    
    std::vector<double> get_Values()
    {
        return solx;
    };
};


void mexFunction (int _no_op_args, mxArray *mex_op[], int _no_ip_args, const mxArray *mex_in[] )
{
    /////////////////// DECLARATIONS ///////////////////////////////////
    // I/O variables
    mxDouble energy[1];
    
    /////////////////// INPUT ///////////////////////////////////
    Eigen::Map<Eigen::VectorXd,Eigen::Aligned> sheet_prop (mxGetPr(mex_in[1]), 5, 1);
    int* no_particles = (int*) mxGetData( mex_in[3] );
    int No_particles = no_particles[0];
    Eigen::Map<Eigen::ArrayXXi,Eigen::Aligned> node_info ((int*)mxGetData(mex_in[2]), no_particles[0], 14);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> coords (mxGetPr(mex_in[0]), no_particles[0], 3);
    Eigen::Map<Eigen::VectorXd,Eigen::Aligned> free_indices (mxGetPr(mex_in[4]), mxGetNumberOfElements(mex_in[4]), 1);
    Eigen::Map<Eigen::VectorXd,Eigen::Aligned> X_init (mxGetPr(mex_in[5]), mxGetNumberOfElements(mex_in[5]), 1);
    
    
    /////////////////// CREATING OPTIMIZATION OBJECT //////////////////
    opt_obj optimizer( coords,sheet_prop,node_info,No_particles,free_indices,X_init );
    bool flag = optimizer.solveOPT;
    Eigen::VectorXd x_opt = optimizer.get_Values();
    
    /////////////////// OUTPUT ///////////////////////////////////
    mex_op[0] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    memcpy(mxGetPr(mex_op[0]), energy, sizeof(double));
    
    return;
}