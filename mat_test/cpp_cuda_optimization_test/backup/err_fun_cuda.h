#ifndef ERROR_FUN_CUDA
#define ERROR_FUN_CUDA

#include <string>
#include <vector>
#include </usr/local/include/eigen3/Eigen/Eigen>

double error_fun_cuda(int, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd, std::vector<double>);

#endif