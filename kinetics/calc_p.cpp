#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void calc_p(double* r, double* A_links_transform, double* p)
{   
    double* a_tempt = new double[3];
    double* r_tempt = new double[3];
    double* A_links_transform_tempt_i = new double[9];
    double* A_links_transform_tempt_i_multi_a_tempt = new double[3];
    
    SetZeroMatrix_( 1,  3,  a_tempt);
    SetZeroMatrix_( 1,  3,  r_tempt);
    SetZeroMatrix_( 1,  9,  A_links_transform_tempt_i);
	
	for(int i = 0; i < N; i++)
	{		
	    SetZeroMatrix_( 1, 3, A_links_transform_tempt_i_multi_a_tempt);
		
		MatrixExtract_( 1, 3*N, 1, 1, i*3+1, i*3+3, a, a_tempt);
		
		ScaleMatrix_(1, 3, (-1),  a_tempt, a_tempt);
		
		MatrixExtract_( 1, 3*N, 1, 1, i*3+1, i*3+3, r, r_tempt);
		
		MatrixExtract_( 1, 3*3*N, 1, 1, i*9+1, i*9+9, A_links_transform, A_links_transform_tempt_i);
		
		/*
		for(int j = 0; j < 3; j++)
		{
			a_tempt[j] = (-1) * a[i * 3 + j];
			r_tempt[j] = r[i * 3 + j];
			
			for(int k = 0; k < 3; k++)
				A_links_transform_tempt_i[j * 3 + k] = A_links_transform[i * 9 + j * 3 + k];
				
		}
		*/

		MatrixMulti_(3, 3, 1, A_links_transform_tempt_i, a_tempt, A_links_transform_tempt_i_multi_a_tempt);
		
		for(int k = 0; k < 3; k++)
			p[i * 3 + k] = r_tempt[k] + A_links_transform_tempt_i_multi_a_tempt[k];
		
	}
    
    delete[] a_tempt;
    delete[] r_tempt;
    delete[] A_links_transform_tempt_i;
    delete[] A_links_transform_tempt_i_multi_a_tempt;
	
}




















