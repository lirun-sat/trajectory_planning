#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

void calc_r(double* r_b, double* A_b, double* A_links_transform, double* r)
{   
    double* a_tempt = new double[3];
    double* b_tempt = new double[3];
    
    SetZeroMatrix_( 1, 3, a_tempt);
    SetZeroMatrix_( 1, 3, b_tempt);
    
    double* A_links_transform_tempt = new double[9];
    double* A_links_transform_tempt_pre = new double[9];
    double* A_links_transform_multi_a_tempt = new double[3];
	double* A_links_transform_tempt_pre_multi_b_tempt = new double[3];
    double* A_b_multi_b_b = new double[3];
    
	
	for(int i = 0; i < N; i++)
	{
        if(i == 0)
        {
            for(int j = 0; j < 3; j++)
                a_tempt[j] = a[i * 3 + j];
            
            for(int j = 0; j< 3; j++)
				for(int k = 0; k < 3; k++)
                    A_links_transform_tempt[j * 3 + k] = A_links_transform[i * 9 + j * 3 + k];
            
            MatrixMulti_(3, 3, 1, A_links_transform_tempt, a_tempt, A_links_transform_multi_a_tempt);
            
            MatrixMulti_(3, 3, 1, A_b, b_b, A_b_multi_b_b);
            
            for(int j = 0; j < 3; j++)
                r[i * 3 + j] = r_b[j] + A_b_multi_b_b[j] + A_links_transform_multi_a_tempt[j];
            
        }
        else
        {
            for(int j = 0; j < 3; j++)
            {
                a_tempt[j] = a[i * 3 + j];
                b_tempt[j] = b[(i-1) * 3 + j];
            }
            
            for(int j = 0; j< 3; j++)
				for(int k = 0; k < 3; k++)
				{
					A_links_transform_tempt[j * 3 + k] = A_links_transform[i * 9 + j * 3 + k];
					A_links_transform_tempt_pre[j * 3 + k] = A_links_transform[(i-1) * 9 + j * 3 + k];
				}
            
            MatrixMulti_(3, 3, 1, A_links_transform_tempt, a_tempt, A_links_transform_multi_a_tempt);
            
            MatrixMulti_(3, 3, 1, A_links_transform_tempt_pre, b_tempt, A_links_transform_tempt_pre_multi_b_tempt);
                        
            for(int j = 0; j < 3; j++)
			    r[i * 3 + j] = r[(i-1) * 3 + j] + A_links_transform_tempt_pre_multi_b_tempt[j] + A_links_transform_multi_a_tempt[j];
            
        }		    
	
	}
    
    delete[] a_tempt;
    delete[] b_tempt;
    delete[] A_links_transform_tempt;
	delete[] A_links_transform_tempt_pre;
	delete[] A_links_transform_multi_a_tempt;
	delete[] A_links_transform_tempt_pre_multi_b_tempt;
	delete[] A_b_multi_b_b;	
	
}


