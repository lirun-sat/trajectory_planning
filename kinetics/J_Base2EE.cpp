
#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void J_Base2EE(double* r_e, double* r_b, double* J_bE)
{
    double* tempt;
    tempt = new double[9];
    
    double* r_b2r_e;
    r_b2r_e = new double[3];
	
	
	for(int i = 0; i < 6; i++)
    {
		for(int j = 0; j < 6; j++)
			{
				if(i == j)
					J_bE[i*6+j] = 1;
		
				else
					J_bE[i*6+j] = 0;
		
			}
    }
	
	// double r_b2r_e[3] = {0};
	
	for(int i = 0; i < 3; i++)
		r_b2r_e[i] = r_e[i] - r_b[i];
	
	cross(r_b2r_e, tempt);
	
	for(int i = 0; i < 9; i++)
		tempt[i] = tempt[i] * (-1);    // 得到 -cross(r_e - r_b)
	
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			J_bE[i * 6 + (3 + j)] = tempt[i * 3 + j];
    
    
    delete[] tempt;
    delete[] r_b2r_e;
		
}
