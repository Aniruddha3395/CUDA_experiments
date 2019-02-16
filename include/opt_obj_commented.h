#ifndef OPT_OBJ_RW
#define OPT_OBJ_RW

#include <nlopt.hpp>
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <vector>

namespace opt_obj
{
    class opt_obj
    {
    private:

    public:
        // constructor and distructor
        opt_obj( Eigen::MatrixXd _coords, Eigen::MatrixXd _prop, Eigen::MatrixXd _nodeInfo, 
            int _NoParticles, Eigen::MatrixXd _freeIndices, std::vector<double> _X_init, double OptH, double OptXtolRel);
        ~opt_obj();

        // other variables
        Eigen::MatrixXd coords;
        Eigen::MatrixXd properties;
        Eigen::MatrixXd node_info;
        int no_particles;
        Eigen::MatrixXd free_indices;
        std::vector<double> X_init;
        // Eigen::MatrixXd spr_pot;
        // Eigen::MatrixXd shr_pot;
        // Eigen::MatrixXd bndg_pot;
        int counter;
        // double *ptr_node_info;
        // double *ptr_coords;
        // double *ptr_spr_pot;
        // double *ptr_shr_pot;
        // double *ptr_bndg_pot;

        int block_count;
        
        // optim solver
        nlopt::opt opt;
        int OptVarDim;
        double optXtolRel;
        double optH;
        std::vector<double> optx;
        std::vector<double> solx;
        double solminf;
        nlopt::algorithm alg_type;
        std::vector<double> OptVarlb;
        std::vector<double> OptVarub;   
        double ErrFun_cuda(const std::vector<double> &x);
        double ObjFun(const std::vector<double> &x, std::vector<double> &grad);
        bool solveOPT();
        double get_solminf();
        std::vector<double> get_solx();
    };
}
#endif