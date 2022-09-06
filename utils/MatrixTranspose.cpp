#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

void MatrixTranspose(double* M)
{
	// double M_tempt[9] = {0};
    double* M_tempt;
    M_tempt = new double[3*3];
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			M_tempt[i * 3 + j] = M[i * 3 + j];
			
	M[0] = M_tempt[0];
	M[1] = M_tempt[3];
	M[2] = M_tempt[6];
	
	M[3] = M_tempt[1];
	M[4] = M_tempt[4];
	M[5] = M_tempt[7];
	
	M[6] = M_tempt[2];
	M[7] = M_tempt[5];
	M[8] = M_tempt[8];
    
    delete[] M_tempt;
	
}
