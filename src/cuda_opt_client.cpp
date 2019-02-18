#include "ros/ros.h"
#include "cuda_opt/OptCuda.h"
#include <cstdlib>

int main(int argc, char **argv)
{
    ros::init(argc, argv, "drape_sim_client");

    ros::NodeHandle n;
    ros::ServiceClient client = n.serviceClient<cuda_opt::OptCuda>("opt_w_cuda");
    // beginner_tutorials::AddTwoInts srv;
    // beginner_tutorials::serv4 srv;
    cuda_opt::OptCuda srv;
    srv.request.file_num = argv[1];
    if (client.call(srv))
    {
        ROS_INFO("error val: %f", (float)srv.response.error);

    }
    else
    {
        ROS_ERROR("Failed to call service opt_w_cuda");
        return 1;
    }
    return 0;
}