#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void MatrixMulti(double* A, double* B, double* AB)
{
	//两个 3x3 矩阵相乘
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			AB[i*3+j] = 0;
	
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


/*
int main()
{
	double xx[9] = {1, 2, 3, 
				   3, 2, 1, 
				   4, 5, 6};
	
	double yy[9] = {1, 1, 1, 
				   1, 1, 1, 
				   1, 1, 1};
	
    double* M = new double[3 * 3];
    MatrixMulti(xx, yy, M);
    
	for(int i = 0; i < (3*3); i++)
		cout << M[i] << endl;
    
    delete[] M;

    return 0;
}

*/



