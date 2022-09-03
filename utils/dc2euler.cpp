
#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

void dc2euler(double* R, double* euler)
{
	double theta_y_1;
	double theta_y_2;
	double theta_x_1;
	double theta_x_2;
	double theta_z_1;
	double theta_z_2;
	
	if(R[2 * 3 + 0] != 1 && R[2 * 3 + 0] != (-1))
	{
		theta_y_1 = -asin(R[2 * 3 + 0]);
		theta_y_2 = M_PI - theta_y_1;
		
		theta_x_1 = atan2(R[2 * 3 + 1] / cos(theta_y_1), R[2 * 3 + 2] / cos(theta_y_1));
		theta_x_2 = atan2(R[2 * 3 + 1] / cos(theta_y_2), R[2 * 3 + 2] / cos(theta_y_2));
		
		theta_z_1 = atan2(R[1 * 3 + 0] / cos(theta_y_1), R[0 * 3 + 0] / cos(theta_y_1));
        theta_z_2 = atan2(R[1 * 3 + 0] / cos(theta_y_2), R[0 * 3 + 0] / cos(theta_y_2));
		
	}
	
	else if(R[2 * 3 + 0] == (-1))
	{
		theta_z_1 = 0;
        theta_z_2 = 0;
		theta_y_1 = M_PI / 2;
		theta_y_2 = M_PI / 2;
		
		theta_x_1 = theta_z_1 + atan2(R[0 * 3 + 1], R[0 * 3 + 2]);
		theta_x_2 = theta_z_1 + atan2(R[0 * 3 + 1], R[0 * 3 + 2]);
		theta_x_1 = theta_z_1 + atan2(R[0 * 3 + 1], R[0 * 3 + 2]);
        theta_x_2 = theta_z_1 + atan2(R[0 * 3 + 1], R[0 * 3 + 2]);
		
	}
	
	else
	{
		theta_z_1 = 0;
        theta_z_2 = 0;
        theta_y_1 = -M_PI / 2;
        theta_y_2 = -M_PI / 2;
        theta_x_1 = -theta_z_1 + atan2(-R[0 * 3 + 1], -R[0 * 3 + 2]);
        theta_x_2 = -theta_z_1 + atan2(-R[0 * 3 + 1], -R[0 * 3 + 2]);
	}
	
    euler[0] = theta_z_1; 
    euler[1] = theta_y_1;
    euler[2] = theta_x_1;
    euler[3] = theta_z_2;
    euler[4] = theta_y_2;
    euler[5] = theta_x_2;
	
}





















