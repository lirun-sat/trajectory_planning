void SetZeroMatrix_(int row, int col, double* A)
{   

    for(int i = 0; i < row; i++)
    	for(int j = 0; j < col; j++)
    		A[i * col + j] = 0;

				
}

