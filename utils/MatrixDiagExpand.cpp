#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

void MatrixDiagExpand(double* A, int row, int col, double* A_expand)
{
	for(int i = 0; i < 2*row; i++)
	{
		for(int j = 0; j < 2*col; j++)
		{
			A_expand[2*col*i + j] = 0;
		}
	}
	
	for(int i = 0; i < 2*row; i++)
	{
		if(i < row)
		{
			for(int j = 0; j < col; j++)
			{
				A_expand[2*col*i + j] = A[col*i + j];
			}
		}
		else
		{
			for(int j = col; j < 2*col; j++)
			{
				A_expand[2*col*i + j] = A[col*(i-row) + (j-col)];
			}
		}
		
	}
	
	
}

/*
int main()
{
	double A[] = {1,2,3,
				  4,5,6};
	
	double *B = new double[2*2 * 2*3];
	
	MatrixDiagExpand(A, 2, 3, B);
	
	for(int i = 0; i < 2*2; i++)
	{
		cout << endl;
		for(int j=0;j<2*3;j++)
		{
			cout << B[2*3*i + j] << "  ";
		}
		cout << endl;
	}
	
}
*/

/*
int main()
{
	double A[] = {1,2,3,
				  4,5,6,
				  7,8,9};
	
	double *B = new double[2*3 * 2*3];
	
	MatrixDiagExpand(A, 3, 3, B);
	
	for(int i = 0; i < 2*3; i++)
	{
		cout << endl;
		for(int j=0;j<2*3;j++)
		{
			cout << B[2*3*i + j] << "  ";
		}
		cout << endl;
	}
	
}
*/

/*
int main()
{
	double A[] = {1,2,
				  4,5,
				  7,8};
	
	double *B = new double[2*3 * 2*2];
	
	MatrixDiagExpand(A, 3, 2, B);
	
	for(int i = 0; i < 2*3; i++)
	{
		cout << endl;
		for(int j=0;j<2*2;j++)
		{
			cout << B[2*2*i + j] << "  ";
		}
		cout << endl;
	}
	
}
*/



