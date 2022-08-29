#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

void MatrixTranspose_(int row, int col, double* M, double* M_T)
{
	// M 是一个 任意 的 m * n 的矩阵, 输出 M 矩阵 的 转置
	
    int m = row;
	int n = col;
    
    double* M_tempt;
    M_tempt = new double[m * n];
    
    double* M_T_tempt;
    M_T_tempt = new double[n * m];
    
	for(int i = 0; i < m; i++)
    	for(int j = 0; j < n; j++)
    		M_tempt[i * n + j] = M[i * n + j];
	
	for(int p = 0; p < n; p++)
		for(int q = 0; q < m; q++)
			M_T_tempt[p * m + q] = M_tempt[q * n + p];
    
    for(int i = 0; i < n; i++)
    	for(int j = 0; j < m; j++)
    		M_T[i * m + j] = M_T_tempt[i * m + j];
    
    delete[] M_tempt;
    delete[] M_T_tempt;
	
}



/*
int main()
{
	double X[] = {1, 2, 3,
				 4, 5, 6,
				 7, 8, 9,
				 1, 4, 7, 
				 3, 9, 5};
	
	double* X_T = new double[15];

    MatrixTranspose_(5, 3, X, X_T);

    cout << X_T[0] << " " << X_T[1] << " " << X_T[2] << " " << X_T[3] << " " << X_T[4]<< endl;
    cout << X_T[5] << " " << X_T[6] << " " << X_T[7] << " " << X_T[8] << " " << X_T[9]<< endl;
	cout << X_T[10] << " " << X_T[11] << " " << X_T[12] << " " << X_T[13] << " " << X_T[14]<< endl;
    
	delete[] X_T;
	
    return 0;
}
*/

