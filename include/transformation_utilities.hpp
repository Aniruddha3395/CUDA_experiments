#ifndef TRANSFORMATION_UTILITIES
#define TRANSFORMATION_UTILITIES

#include <Eigen/Eigen>
#include <vector>
#include <string>

class rtf
{
	// NOTE: For inputs nad outputs-
	// Rotation Matrix [3x3]
	// Euler Angles [1x3]
	// Quaternion [1x4]

public:

	// Make homogenous transformation matrix from translation vector and rotation matrix 
	static Eigen::Matrix4i hom_T(Eigen::Vector3i, Eigen::Matrix3i);
	static Eigen::Matrix4f hom_T(Eigen::Vector3f, Eigen::Matrix3f);
	static Eigen::Matrix4d hom_T(Eigen::Vector3d, Eigen::Matrix3d);

	// apply homogenous transformation to data points
	static Eigen::MatrixXd apply_transformation(Eigen::MatrixXd, Eigen::Matrix4d);

	// Verifies validity of a rotation sequence
    static std::string validate_seq(std::string="");
	
	// Euler Angles to Rotation Matrix Conversion 
	// NOTE: Euler Angles are: (alpha,beta,gamma), Default sequence: ZYX 
	static Eigen::Matrix3d eul2rot(Eigen::MatrixXd, std::string="");

	// Rotation Matrix to Euler Angles Conversion 
	// NOTE: Euler Angles are: (alpha,beta,gamma), Default sequence: ZYX
	static Eigen::MatrixXd rot2eul(Eigen::Matrix3d, std::string="");

	// Quaternion to Rotation Matrix Conversion
	// NOTE: Input Quaternion is in [x,y,z,w] form  
	static Eigen::Matrix3d qt2rot(Eigen::MatrixXd);

	// Rotation Matrix to Quaternion Conversion
	// NOTE: Input Quaternion is in [x,y,z,w] form 
	static Eigen::MatrixXd rot2qt(Eigen::Matrix3d);

	// Euler Angles to Quaternion Matrix Conversion 
	// NOTE: Euler Angles are: (alpha,beta,gamma), Default sequence: ZYX 
	// NOTE: Input Quaternion is in [x,y,z,w] form 
	static Eigen::MatrixXd eul2qt(Eigen::MatrixXd, std::string="");

	// Quaternion to Euler Angles Matrix Conversion 
	// NOTE: Euler Angles are: (alpha,beta,gamma), Default sequence: ZYX
	// NOTE: Input Quaternion is in [x,y,z,w] form 
	static Eigen::MatrixXd qt2eul(Eigen::MatrixXd, std::string="");

	// Euler Angles to Rotation Vectors Conversion
	// Input euler alpha,beta,gamma for ZYX (in radians)
	static Eigen::MatrixXd eul2bxbybz(Eigen::MatrixXd);

	// Euler Angles to Rotation Vectors Conversion
	// Input bxbybz as a single nx9 Matrix
	// Output is ZYX euler angles in radians 	
	static Eigen::MatrixXd bxbybz2eul(Eigen::MatrixXd);

	// Generate homogenous transformation between 2 coordiante frames using known correspondances between 4 or more points 
	static Eigen::MatrixXd get_rob_T_part(Eigen::MatrixXd, Eigen::MatrixXd);
	static Eigen::MatrixXd mean(Eigen::MatrixXi);
	static Eigen::MatrixXd mean(Eigen::MatrixXf);
	static Eigen::MatrixXd mean(Eigen::MatrixXd);
	
	// Generate rotation matrix about X-axis
	// NOTE: input angle must be radians
	static Eigen::Matrix3d rot_x(double);

	// Generate rotation matrix about Y-axis
	// NOTE: input angle must be radians
	static Eigen::Matrix3d rot_y(double);

	// Generate rotation matrix about Z-axis
	// NOTE: input angle must be radians
	static Eigen::Matrix3d rot_z(double);
};

#endif