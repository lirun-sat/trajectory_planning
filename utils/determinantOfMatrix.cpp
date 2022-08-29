#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

double calc_determinantOfMatrix(double* a, int num, int row, int col, int n)
{
    /*
     n 表示 要求行列式的矩阵
     a 是一个一维数组；
     num 是数组 a 的元素个数；
     row 是 a 要转换成的二维数组的 行数；
     col 是 a 要转换成的二维数组的 列数；
     该程序中，因为是求 矩阵的行列式 ， 所以 row = col = n, num = n**2
    */
    double num1 =1;
	double num2 = 1;
	double det = 1;
	int index = 1;
	double total = 1; // Initialize result
    double result_tempt = 0;
    double result = 0;
        
    double** mat_tempt;
    mat_tempt = new double* [n];
    
    // temporary array for storing row
    double* temp;
    temp = new double[n];    
    
    for(int k = 0; k < row; k++)
    {
        mat_tempt[k] = new double[col];
        
        for(int j = 0; j < (num / row); j++)
            mat_tempt[k][j] = a[k * col + j];
    }
    
    
	// loop for traversing the diagonal elements
	for (int i = 0; i < n; i++)
	{
		index = i; // initialize the index

		// finding the index which has non zero value
		while (index < n && mat_tempt[index][i] == 0)
		{
			index++;
		}
		if (index == n) // if there is non zero element
		{
			// the determinant of matrix as zero
			continue;
		}
		if (index != i)
		{
			// loop for swapping the diagonal element row and
			// index row
			for (int j = 0; j < n; j++)
			{
				swap(mat_tempt[index][j], mat_tempt[i][j]);
			}
			// determinant sign changes when we shift rows
			// go through determinant properties
			det = det * pow(-1, index - i);
		}

		// storing the values of diagonal row elements
		for (int j = 0; j < n; j++)
		{
			temp[j] = mat_tempt[i][j];
		}
		// traversing every row below the diagonal element
		for (int j = i + 1; j < n; j++)
		{
			num1 = temp[i]; // value of diagonal element
			num2 = mat_tempt[j][i]; // value of next row element

			// traversing every column of row
			// and multiplying to every row
			for (int k = 0; k < n; k++)
			{
				// multiplying to make the diagonal
				// element and next row element equal
				mat_tempt[j][k] = (num1 * mat_tempt[j][k]) - (num2 * temp[k]);
                
			}
            
			total = total * num1; // Det(kA)=kDet(A);
            
		}
	}

	// multiplying the diagonal elements to get determinant
	for (int i = 0; i < n; i++)
		det = det * mat_tempt[i][i];
    
    
    result_tempt = det / total;
    result = result_tempt;
    
    
    for (int i = 0; i < n; i++)
        delete[] mat_tempt[i];
    
    delete[] mat_tempt;
    
    delete[] temp;
    
    return result;
		
}


/*
// Driver code
int main()
{
	//double mat_tempt[10][10];
	
	//double a = 1;
	
	//double mat[10][10];
	//for (int i = 0; i < 2; i++)
    //    for (int j = 0; j < 2; j++)
     //       mat[i][j] = a++;
	
	// double* a;
	double a[16] = {1, 3, 5, 7,
			 9, 3, 4, 56,
		     6, 32, 9, 6,
		     5, 2, 1, 8 
		  };
					    
			
	int num = sizeof(a) / sizeof(a[0]);
    printf("num is : %d \n", num);
    
    int n = sqrt (num);
    
    double result;
    
    result = calc_determinantOfMatrix(a, num, n,  n,  n);

	// Function call
	printf("Determinant of the matrix is : %f \n", result);

		
	return 0;
}

*/




