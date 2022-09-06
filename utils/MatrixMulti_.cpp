#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

void MatrixMulti_(int row_1, int col_1, int col_2, double* A, double* B, double* AB)
{
	
	/*
	 两个 任意维数 矩阵相乘
	 row_1 第一个矩阵的行数
	 col_1 第一个矩阵的列数
	 col_1 第二个矩阵的行数和第一个矩阵的列数相等
	 col_2 第二个矩阵的列数
	*/
    
    double* A_tempt;
    A_tempt = new double[row_1 * col_1];
    
    double* B_tempt;
    B_tempt = new double[col_1 * col_2];
    
    double* AB_tempt;
    AB_tempt = new double[row_1 * col_2];
    
    
    for(int i = 0; i < row_1; i++)
    	for(int j = 0; j < col_1; j++)
    		A_tempt[i * col_1 + j] = A[i * col_1 + j];
    
    for(int i = 0; i < col_1; i++)
    	for(int j = 0; j < col_2; j++)
    		B_tempt[i * col_2 + j] = B[i * col_2 + j];
    
    
    for(int i = 0; i < row_1; i++)
    	for(int j = 0; j < col_2; j++)
    		AB_tempt[i * col_2 + j] = 0;
    		
    
	for(int i = 0; i < row_1; i++)
		for(int j = 0; j < col_2; j++)
			for(int k = 0; k < col_1; k++)
				AB_tempt[i * col_2 + j] += A_tempt[i * col_1 + k] * B_tempt[k * col_2 + j];
    
    for(int i = 0; i < row_1; i++)
    	for(int j = 0; j < col_2; j++)
    		AB[i * col_2 + j] = AB_tempt[i * col_2 + j];
    
    delete[] A_tempt;
    delete[] B_tempt;
    delete[] AB_tempt;
    
				
}


/*
int main()
{
	double xx[] = {1, 2, 3, 
				  3, 4, 5,
				  5, 6, 2};
	
	double yy[] = {4, 3, 5};
	
    double* M = new double[3];
	
    MatrixMulti_(3, 3, 1, xx, yy, M);
    
	for(int i = 0; i < (3); i++)
		cout << M[i] << endl;
    
    delete[] M;

    return 0;
}


*/
