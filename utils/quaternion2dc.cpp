#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void quaternion2dc(double q_0, double q_1, double q_2, double q_3, double* dc)
{	
	dc[0] = pow(q_0, 2) + pow(q_1, 2) - pow(q_2, 2) - pow(q_3, 2);
	dc[1] = 2 * (q_1 * q_2 - q_0 * q_3);
	dc[2] = 2 * (q_1 * q_3 + q_0 * q_2);
	
	dc[3] = 2 * (q_1 * q_2 + q_0 * q_3);
	dc[4] = pow(q_0, 2) - pow(q_1, 2) + pow(q_2, 2) - pow(q_3, 2);
	dc[5] = 2 * (q_2 * q_3 - q_0 * q_1);
	
	dc[6] = 2 * (q_1 * q_3 - q_0 * q_2);
	dc[7] = 2 * (q_2 * q_3 + q_0 * q_1);
	dc[8] = pow(q_0, 2) - pow(q_1, 2) - pow(q_2, 2) + pow(q_3, 2);
    
}
