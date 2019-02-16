#ifndef UTILITIES
#define UTILITIES

#include <Eigen/Dense>
#include <vector>
#include <string>

class ut
{
public:
	// get unique rows in a vector
	static std::vector<std::vector<int> > GetUniqueRows(std::vector<std::vector<int> >);
	static std::vector<std::vector<float> > GetUniqueRows(std::vector<std::vector<float> >);
	static std::vector<std::vector<double> > GetUniqueRows(std::vector<std::vector<double> >);
	
	// vector to matrix conversion
	static Eigen::MatrixXi vec_to_mat(std::vector<std::vector<int> >);
	static Eigen::MatrixXf vec_to_mat(std::vector<std::vector<float> >);
	static Eigen::MatrixXd vec_to_mat(std::vector<std::vector<double> >);
	
	
	// matrix to vector conversion
	static std::vector<std::vector<int> > mat_to_vec(Eigen::MatrixXi);
	static std::vector<std::vector<float> > mat_to_vec(Eigen::MatrixXf);
	static std::vector<std::vector<double> > mat_to_vec(Eigen::MatrixXd);

	// print vector
	static void disp_vec(std::vector<std::vector<int> >);
	static void disp_vec(std::vector<std::vector<float> >);
	static void disp_vec(std::vector<std::vector<double> >);
	
	// Compute bxbybz (nx9 vector) from normals
	static Eigen::MatrixXd compute_TCP(Eigen::MatrixXd, Eigen::MatrixXd);
	static double get_pt_to_lsf_plane_dist(Eigen::MatrixXd, Eigen::MatrixXd);
	static Eigen::MatrixXd get_traj_wrt_tcp(Eigen::Matrix4d, std::vector<std::vector<double> >);

	// sort matrix based on the input column
	static std::vector<std::vector<int> > SortRows(std::vector<std::vector<int> >, int);
	static std::vector<std::vector<float> > SortRows(std::vector<std::vector<float> >, int);
	static std::vector<std::vector<double> > SortRows(std::vector<std::vector<double> >, int);
	
	// check if row/vector is a member of matrix/2D vector
	static std::vector<int> ismember(std::vector<std::vector<int> >, std::vector<int>);
	static std::vector<int> ismember(std::vector<std::vector<float> >, std::vector<float>);
	static std::vector<int> ismember(std::vector<std::vector<double> >, std::vector<double>);
	static std::vector<int> ismember(Eigen::MatrixXi, Eigen::MatrixXi);
	static std::vector<int> ismember(Eigen::MatrixXf, Eigen::MatrixXf);
	static std::vector<int> ismember(Eigen::MatrixXd, Eigen::MatrixXd);

	// find rowwise mean
	static Eigen::MatrixXd mean(Eigen::MatrixXi);
	static Eigen::MatrixXd mean(Eigen::MatrixXf);
	static Eigen::MatrixXd mean(Eigen::MatrixXd);
	
	// find median
	static double median(std::vector<int>);
	static double median(std::vector<float>);
	static double median(std::vector<double>);

	// linearly spaced points between start and end point
	static Eigen::VectorXd linsp(double, double, double);

	// inpolygon function to check if point/s is/are inside polygon 
	static std::vector<int> InPoly(Eigen::MatrixXd, Eigen::MatrixXd);

	// find the index of non-zero elements
	static std::vector<int> find_idx(Eigen::VectorXi);
	static std::vector<int> find_idx(Eigen::VectorXf);
	static std::vector<int> find_idx(Eigen::VectorXd);
	static std::vector<int> find_idx(std::vector<int>);
	static std::vector<int> find_idx(std::vector<float>);
	static std::vector<int> find_idx(std::vector<double>);
};

#endif
