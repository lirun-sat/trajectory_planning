void ScaleMatrix_(int row, int col, double scale, double* A, double* scale_A)
{   
    double* A_tempt;
    A_tempt = new double[row * col];
    
    double* AB_tempt;
    AB_tempt = new double[row * col];
    
    
    for(int i = 0; i < row; i++)
    	for(int j = 0; j < col; j++)
    		A_tempt[i * col + j] = A[i * col + j];


    for(int i = 0; i < row; i++)
    	for(int j = 0; j < col; j++)
    		AB_tempt[i * col + j] = 0;
    		
    
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			AB_tempt[i * col + j] = scale * A_tempt[i * col + j];
    
    
    for(int i = 0; i < row; i++)
    	for(int j = 0; j < col; j++)
    		scale_A[i * col + j] = AB_tempt[i * col + j];
    
    delete[] A_tempt;
    delete[] AB_tempt;
    
				
}

