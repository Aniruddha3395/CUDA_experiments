#include <iostream>
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <vector>
#include <string>
#include <cmath>
#include <nlopt.hpp>
#include "stdlib.h"
#include "opt_obj.hpp"

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
    opt_obj::opt_obj( Eigen::MatrixXd _coords, Eigen::MatrixXd _prop, Eigen::MatrixXd _nodeInfo,
                int _NoParticles, Eigen::MatrixXd _freeIndices, std::vector<double> _X_init, double OptH, double OptXtolRel) 
    {

        //choose optimizer
        // optalg = nlopt::LN_NEWUOA;
        // optalg = nlopt::LN_NEWUOA_BOUND;
        // optalg = nlopt::LN_BOBYQA;
        // optalg = nlopt::LN_COBYLA;
        alg_type = nlopt::LD_SLSQP;
        // alg_type = nlopt::LD_LBFGS;
        // optalg = nlopt::GN_ISRES;
        
        
        coords = _coords;
        properties = _prop;
        node_info = _nodeInfo;
        no_particles = _NoParticles;
        free_indices = _freeIndices;
        X_init = _X_init;
        optH = OptH;
        optXtolRel = OptXtolRel;
        solminf = 0;

        // Optimization Parameters
        OptVarDim = _X_init.size();
        opt = nlopt::opt(alg_type, OptVarDim);
        OptVarlb.resize(OptVarDim);
        OptVarub.resize(OptVarDim);
        
        for (int i=0;i<OptVarlb.size();++i)
        {
            OptVarlb[i] = -100000;
            OptVarub[i] = 100000;
        }

        opt.set_xtol_rel(optXtolRel);
        opt.set_min_objective(customminfunc, this);

        // opt.set_lower_bounds(OptVarlb);
        // opt.set_upper_bounds(OptVarub);
        opt.set_maxeval(100);
        

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

    double opt_obj::ErrFun(const std::vector<double> &x)
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

        
        for ( int i = 0; i < free_indices.rows(); i+=3)
        {
            coords(free_indices(i,0)-1,0) = x[i];
            coords(free_indices(i,0)-1,1) = x[i+1];
            coords(free_indices(i,0)-1,2) = x[i+2];
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

    bool opt_obj::solveOPT()
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
        double err = ErrFun(x);
        if (!grad.empty()) {
            std::vector<double> xph = x;
            for (uint i=0; i < x.size(); ++i)
            {
                xph[i] += optH;
                grad[i] = (ErrFun(xph)-err)/optH;
                xph[i] -= optH;
            }
        }    
        return err;
    };
}