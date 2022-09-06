
#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void links_transform(double* A_base, double* q, double* A_links_transform)
{
	
	double* tempt_rpy2dc_joints_initial = new double[3*3];
	double* tempt_rpy2dc_joints_var = new double[3*3];
	double* tempt_A_links_transform_pre = new double[3*3];
	double* tempt_A_links_transform_curr = new double[3*3];
	double* tempt = new double[3*3];
	
	for(int i = 0; i < N; i++)
	{
		if (i == 0)
		{
			rpy2dc(rpy_joints[i * 3], rpy_joints[i * 3 + 1], rpy_joints[i * 3 + 2], tempt_rpy2dc_joints_initial);

			MatrixMulti_(3, 3, 3, A_base, tempt_rpy2dc_joints_initial, tempt);
			
			rpy2dc(0, 0, q[i], tempt_rpy2dc_joints_var);
			
			MatrixMulti_(3, 3, 3, tempt, tempt_rpy2dc_joints_var, tempt_A_links_transform_curr);
			
			for(int j = 0; j< 3; j++)
            {
				for(int k = 0; k < 3; k++)
				{
					A_links_transform[i * 9 + j * 3 + k] = tempt_A_links_transform_curr[j * 3 + k];

				}
            }

		}
		
		else
		{
			rpy2dc(rpy_joints[i * 3], rpy_joints[i * 3 + 1], rpy_joints[i * 3 + 2], tempt_rpy2dc_joints_initial);
			
			for(int j = 0; j < 3; j++)
				for(int k = 0; k < 3; k++)
					tempt_A_links_transform_pre[j * 3 + k] = A_links_transform[(i-1) * 9 + j * 3 + k];
			
			MatrixMulti_(3, 3, 3, tempt_A_links_transform_pre, tempt_rpy2dc_joints_initial, tempt);
			
			rpy2dc(0, 0, q[i], tempt_rpy2dc_joints_var);
			
			MatrixMulti_(3, 3, 3, tempt, tempt_rpy2dc_joints_var, tempt_A_links_transform_curr);
	
			for(int j = 0; j< 3; j++)
            {
				for(int k = 0; k < 3; k++)
				{
					A_links_transform[i * 9 + j * 3 + k] = tempt_A_links_transform_curr[j * 3 + k];

				}
            }

		}
		
	}
	
	delete[] tempt_rpy2dc_joints_initial;
	delete[] tempt_rpy2dc_joints_var;
	delete[] tempt_A_links_transform_pre;
	delete[] tempt_A_links_transform_curr;	
	delete[] tempt;
    
}





/*

int main()
{
    int N = 7;	
	
	double A_base[9] = {1, 0, 0, 
					  0, 1, 0, 
					  0, 0, 1};
	
	double rpy_joints[7*3] = {   0,    0, PI/2,
					        -PI/2, PI/2,    0,
					         PI/2,    0,   PI,
					        -PI/2,    0,    0,
					         PI/2,    0,   PI,
					         PI/2,    0,    0,
					        -PI/2,    0,    0};
	
	// double q[N] = {0, 0, 0, 0, 0, 0, 0};
	double q[N] = {0.03384647, 0.26985292, 0.5218816, 1.87610337, 0.40879103, 1.04341179, -0.80001746};
	
    double* M = new double[N * 3 * 3];
	
    links_transform(A_base, rpy_joints, q, N, M);
    
	for(int i = 0; i < (N*3*3); i++)
		cout << M[i] << endl;
    
    delete[] M;

    return 0;
}
*/
