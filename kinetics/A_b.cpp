
#define _USE_MATH_DEFINES 
#include"utils.h"


void A_b(double alpha_base, double beta_base, double gamma_base, double* A_base)
{
    rpy2dc(alpha_base, beta_base, gamma_base, A_base);
    
}

/*

void rpy2dc(double roll, double pitch, double yaw, double* direction_cosine)
{
    double* tempt_yaw = new double[9];
	double* tempt_pitch = new double[9];
	double* tempt_roll = new double[9];
	double* tempt_yawMulpitch = new double[9];
	
	cxyz(yaw, 0, 0, 1, tempt_yaw);
	cxyz(pitch, 0, 1, 0, tempt_pitch);
	cxyz(roll, 1, 0, 0, tempt_roll);
	
	MatrixMulti(tempt_yaw, tempt_pitch, tempt_yawMulpitch);
	
	MatrixMulti(tempt_yawMulpitch, tempt_roll, direction_cosine);
	
	delete[] tempt_yaw;
	delete[] tempt_pitch;
	delete[] tempt_roll;
	delete[] tempt_yawMulpitch;
    
}

void MatrixMulti(double* A, double* B, double* AB)
{
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if (i == 0)
				AB[i + j] = A[i] * B[j] + A[i + 1] * B[j + 3] + A[i + 2] * B[j + 6];
			
			else if (i == 1)
				AB[i + 2 + j] = A[i + 2] * B[j] + A[i + 2 + 1] * B[j + 3] + A[i + 2 + 2] * B[j + 6];
			
			else
				AB[i + 2 + 2 + j] = A[i + 2 + 2] * B[j] + A[i + 2 + 2 + 1] * B[j + 3] + A[i + 2 + 2 + 2] * B[j + 6];
			
		}
		
	}
				
}

void cxyz(double theta, int x, int y, int z, double* direction_cosines)
{
    if (x == 1)
	{
		direction_cosines[0] = 1;
		direction_cosines[1] = 0;
		direction_cosines[2] = 0;
		
		direction_cosines[3] = 0;
		direction_cosines[4] = cos(theta);
		direction_cosines[5] = -sin(theta);
		
		direction_cosines[6] = 0;
		direction_cosines[7] = sin(theta);
		direction_cosines[8] = cos(theta);
		
	}
	else if (y == 1)
	{
		direction_cosines[0] = cos(theta);
		direction_cosines[1] = 0;
		direction_cosines[2] = sin(theta);
		
		direction_cosines[3] = 0;
		direction_cosines[4] = 1;
		direction_cosines[5] = 0;
		
		direction_cosines[6] = -sin(theta);
		direction_cosines[7] = 0;
		direction_cosines[8] = cos(theta);
		
	}
	else if (z == 1)
	{
		direction_cosines[0] = cos(theta);
		direction_cosines[1] = -sin(theta);
		direction_cosines[2] = 0;
		
		direction_cosines[3] = sin(theta);
		direction_cosines[4] = cos(theta);
		direction_cosines[5] = 0;
		
		direction_cosines[6] = 0;
		direction_cosines[7] = 0;
		direction_cosines[8] = 1;
		
	}
    
}


int main()
{
	double roll = 30 * 3.14 / 180;
	double pitch = 30 * 3.14 / 180;
	double yaw = 30 * 3.14 / 180;
	
    double* M = new double[3 * 3];
    A_b(roll, pitch, yaw, M);
    
	for(int i = 0; i < (3*3); i++)
		cout << M[i] << endl;
    
    delete[] M;

    return 0;
}

*/
