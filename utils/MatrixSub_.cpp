void MatrixSub_(int row, int col, double* A, double* B, double* AB)
{   
    double* A_tempt;
    A_tempt = new double[row * col];
    
    double* B_tempt;
    B_tempt = new double[row * col];
    
    double* AB_tempt;
    AB_tempt = new double[row * col];
    
    
    for(int i = 0; i < row; i++)
    	for(int j = 0; j < col; j++)
    		A_tempt[i * col + j] = A[i * col + j];
    
    for(int i = 0; i < row; i++)
    	for(int j = 0; j < col; j++)
    		B_tempt[i * col + j] = B[i * col + j];
    
    
    for(int i = 0; i < row; i++)
    	for(int j = 0; j < col; j++)
    		AB_tempt[i * col + j] = 0;
    		
    
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			AB_tempt[i * col + j] = A_tempt[i * col + j] - B_tempt[i * col + j];
    
    for(int i = 0; i < row; i++)
    	for(int j = 0; j < col; j++)
    		AB[i * col + j] = AB_tempt[i * col + j];
    
    delete[] A_tempt;
    delete[] B_tempt;
    delete[] AB_tempt;
    
				
}

