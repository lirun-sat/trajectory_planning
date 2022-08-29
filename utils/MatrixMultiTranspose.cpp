#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void MatrixMultiTranspose(double* A, double* B, double* AB)
{
	// 矩阵 A 乘以 B.T， 两个矩阵都是 3 * 3 
	
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			AB[i * 3 + j] = 0;
	
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			AB[i * 3 + j] = A[i * 3] * B[j * 3] + A[i * 3 + 1] * B[j * 3 + 1] + A[i * 3 + 2] * B[j * 3 + 2];
		}
		
	}
				
}


/*
int main()
{
	double xx[9] = {1, 2, 3, 
				   3, 2, 1, 
				   4, 5, 6};
	
	double yy[9] = {2, 1, 2, 
				   1, 1, 1, 
				   1, 2, 2};
	
    double* M = new double[3 * 3];
    MatrixMultiTranspose(xx, yy, M);
    
	for(int i = 0; i < (3*3); i++)
		cout << M[i] << endl;
    
    delete[] M;

    return 0;
}

*/



