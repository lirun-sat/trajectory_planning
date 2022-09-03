
#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void forward_kin(double* para, double* eta_end, double* xi_end, double* Pe, double* eta_b, double* xi_b)
{
    double* q = new double[N];
    double* q_dot = new double[N];

    double* quaternion_base = new double[4];
    double* A_b = new double[3*3];
    double* xi_b_tempt = new double[3];

    double* xi_end_tempt = new double[3];
    double* A_links_transform  = new double[3*3*N];
    double* A_b_multi_b_b = new double[3];
    double* r_b_tempt_1 = new double[3];
    double* r_b_tempt_2 = new double[3];
    double* r_b_tempt_2_temp = new double[3];
    double* a_tempt_i = new double[3];
    double* r_b_tempt_3 = new double[3];
    double* r_b = new double[3];
    double* a_tempt_k = new double[3];
    double* b_tempt_k = new double[3];
    double* l_tempt_k = new double[3];    
    double* A_links_transform_tempt_i = new double[3*3];
    double* A_links_transform_tempt_i_multi_a_tempt_i = new double[3];
    double* A_links_transform_tempt_k = new double[3*3];
    double* A_links_transform_tempt_k_multi_l_tempt_k = new double[3];
    double* r_b_tempt_3_tempt = new double[3];
    double* r = new double[3 * N];
    double* r_e = new double[3];
    double* Pe_initial = new double[3];
    double* A_links_transform_N_minus_1 = new double[3*3];
    double* r_e_tempt = new double[3];
    double* b_tempt_N_minus_1 = new double[3];
    double* p = new double[N*3];
    double* A_links_transform_tempt_i_multi_Ez = new double[3];
    double* A_links_transform_tempt_i_multi_Ez_cross = new double[9];
    double* Jm_v = new double[3*N];
    double* Jm_v_tempt = new double[N*3];
    double* Jm_v_tempt_i = new double[3];
    double* Jm_w = new double[3*N];
    double* Jm_w_tempt = new double[N*3];
    double* p_tempt = new double[3];
    double* r_e_minus_p_i = new double[3];
    double* Jm = new double[6*N];
    double* Jm_tempt = new double[N*6];
    double* J_bm_w = new double[3*N];
    double* J_bm_v = new double[3*N];
    double* J_bm = new double[6*N];	
    double* J_bE = new double[6*6];
    double* J_g = new double[6*N];
    double* J_bE_multi_J_bm = new double[6*N];
    double* J_g_v = new double[3*N];
    double* J_g_w = new double[3*N];    
    double* quaternion_end  = new double[4];
    double eta_b_dot;
    double* eta_b_dot_tempt_1 = new double[N];
    double* eta_b_dot_tempt_2 = new double;  
    double* xi_b_dot = new double[3];
    double* cross_xi_b = new double[3*3];  
    double* xi_b_dot_tempt_1 = new double[3*3];
    double* xi_b_dot_tempt_2 = new double[3*N];
    double* xi_b_dot_tempt_3 = new double[3]; 
    double eta_end_dot;
    double* eta_end_dot_tempt_1 = new double[N];
    double* eta_end_dot_tempt_2 = new double;
    double* xi_end_dot = new double[3];
    double* cross_xi_end = new double[3*3]; 
    double* xi_end_dot_tempt_1 = new double[3*3];
    double* xi_end_dot_tempt_2 = new double[3*N];
    double* xi_end_dot_tempt_3 = new double[3];
    double* v_e = new double[3];
    double* v_b = new double[3];
    double total_mass_for_links = 0;
	double total_mass = 0;	
	double q_dot_temp;
	
    
    for(int i = 0; i < N; i++)
        q[i] = q_INITIAL[i];
    
    zyx2quaternion(RPY_BASE_INITIAL[2], RPY_BASE_INITIAL[1], RPY_BASE_INITIAL[0], quaternion_base);
    quaternion2dc(quaternion_base[0], quaternion_base[1], quaternion_base[2], quaternion_base[3], A_b);
    *eta_b = quaternion_base[0];
    
    for(int i=0;i<3;i++)
        xi_b[i] = quaternion_base[i+1];
    
// ********************************************************************************************************************************************************************************        
    
    
    links_transform(A_b, q, A_links_transform);
      
    for(int i=0;i<3;i++)
        r_b_tempt_2[i] = 0;
    
    for(int i=0;i<3;i++)
        r_b_tempt_3[i] = 0;
    
    // 自由漂浮，初始动量为 0 ， 以 系统质心为 惯性系原点。
    for (int j = 0; j < N; j++)
    {
        total_mass_for_links += m[j];
    }
	total_mass = m_b + total_mass_for_links;
	MatrixMulti_(3, 3, 1, A_b, b_b, A_b_multi_b_b);	
	ScaleMatrix_( 3, 1, total_mass_for_links, A_b_multi_b_b, r_b_tempt_1);
	
	for(int i=0;i<N;i++)
	{
		
		for(int j=0;j<3;j++)
		{
			a_tempt_i[j] = a[i*3+j];			
			r_b_tempt_3_tempt[j] = 0;			
            for (int k = 0; k < 3; k++)
            {
                A_links_transform_tempt_i[j * 3 + k] = A_links_transform[i * 9 + j * 3 + k];
            }
		}
		
		MatrixMulti_(3, 3, 1, A_links_transform_tempt_i, a_tempt_i, A_links_transform_tempt_i_multi_a_tempt_i);		
		ScaleMatrix_( 1, 3, m[i], A_links_transform_tempt_i_multi_a_tempt_i, r_b_tempt_2_temp);		
		MatrixAdd_( 1, 3, r_b_tempt_2, r_b_tempt_2_temp, r_b_tempt_2);

		for(int k=0;k<i;k++)
		{
			for(int jj=0;jj<3;jj++)
			{
				a_tempt_k[jj] = a[k*3+jj];
				b_tempt_k[jj] = b[k*3+jj];
				l_tempt_k[jj] = a_tempt_k[jj] + b_tempt_k[jj];				
                for (int kk = 0; kk < 3; kk++)
                {
                    A_links_transform_tempt_k[jj * 3 + kk] = A_links_transform[k * 9 + jj * 3 + kk];
                }
			}						
			MatrixMulti_(3, 3, 1, A_links_transform_tempt_k, l_tempt_k, A_links_transform_tempt_k_multi_l_tempt_k);			
			MatrixAdd_( 1, 3, r_b_tempt_3_tempt, A_links_transform_tempt_k_multi_l_tempt_k, r_b_tempt_3_tempt);	
		}
		
		ScaleMatrix_( 1, 3, m[i], r_b_tempt_3_tempt, r_b_tempt_3_tempt);		
		MatrixAdd_( 1, 3, r_b_tempt_3, r_b_tempt_3_tempt, r_b_tempt_3);
	}
    //  计算 基座 位置 
    MatrixAdd_( 1, 3, r_b_tempt_1, r_b_tempt_2, r_b_tempt_2);
	MatrixAdd_( 1, 3, r_b_tempt_2, r_b_tempt_3, r_b_tempt_3);	
	ScaleMatrix_( 1, 3, (-1) / total_mass, r_b_tempt_3, r_b);	
	// 计算各个连杆位置	
	calc_r(r_b, A_b, A_links_transform, r);
	for(int i=0;i<3;i++)
	{
		b_tempt_N_minus_1[i] = b[(N-1)*3+i];		
        for (int j = 0; j < 3; j++)
        {
            A_links_transform_N_minus_1[i * 3 + j] = A_links_transform[(N - 1) * 9 + i * 3 + j];
        }
	}
	MatrixMulti_(3, 3, 1, A_links_transform_N_minus_1, b_tempt_N_minus_1, r_e_tempt);
    //  计算末端位置
    for (int i = 0; i < 3; i++)
    {
        r_e[i] = r[(N - 1) * 3 + i] + r_e_tempt[i];
    }
// ********************************************** 初始末端位姿 ************************************************************************************************************   	
    for (int i = 0; i < 3; i++)
    {
        Pe_initial[i] = r_e[i];
    }
    zyx2quaternion(RPY_END_INITIAL[2], RPY_END_INITIAL[1], RPY_END_INITIAL[0], quaternion_end);  
    *eta_end = quaternion_end[0];
    for (int i = 0; i < 3; i++)
    {
        xi_end[i] = quaternion_end[i + 1];
    }
    for (int i = 0; i < 3; i++)
    {
        Pe[i] = Pe_initial[i];
    }
        
// ********************************************************************************************************************************************************************************         
    // 计算 关节位置
	calc_p(r, A_links_transform, p);    

// *************************************  计算 Jm_v, Jm_w, Jm, J_bm_w, J_bm_v, J_bm, J_g_v, J_g_w, J_g  **************************************************************************
    for(int i=0;i<N;i++)
	{

		MatrixExtract_( 1, 3*N, 1, 1, i*3+1, i*3+3, p, p_tempt);		
		ScaleMatrix_( 1, 3, (-1), p_tempt, p_tempt);		
		MatrixAdd_( 1, 3, r_e, p_tempt, r_e_minus_p_i);		
		MatrixExtract_( 1, 3*3*N, 1, 1, i*9+1, i*9+9, A_links_transform, A_links_transform_tempt_i);						
		MatrixMulti_(3, 3, 1, A_links_transform_tempt_i, Ez, A_links_transform_tempt_i_multi_Ez);		
		cross(A_links_transform_tempt_i_multi_Ez, A_links_transform_tempt_i_multi_Ez_cross);		
		MatrixMulti_(3, 3, 1, A_links_transform_tempt_i_multi_Ez_cross, r_e_minus_p_i, Jm_v_tempt_i);		
		
		for(int k = 0; k < 3; k++)
		{
			Jm_v_tempt[i * 3 + k] = Jm_v_tempt_i[k];
			Jm_w_tempt[i * 3 + k] = A_links_transform_tempt_i_multi_Ez[k];
		}
		
		for(int jj = 0; jj < 6; jj++)
		{
			if(jj<3)
				Jm_tempt[i * 6 + jj] = Jm_v_tempt[i * 3 + jj];
			else
				Jm_tempt[i * 6 + jj] = Jm_w_tempt[i * 3 + (jj - 3)];
		}
	}
	MatrixTranspose_(N, 3, Jm_v_tempt, Jm_v);
	MatrixTranspose_(N, 3, Jm_w_tempt, Jm_w);
	MatrixTranspose_(N, 6, Jm_tempt, Jm);

    calc_J_bm( r, r_b, A_b, A_links_transform, p, J_bm_w, J_bm_v);

	for(int j = 0; j < 6; j++)
    {
		for(int k = 0; k < N; k++)
			{
				if(j < 3)
					J_bm[j * N + k] = J_bm_v[j * N + k];
				else
					J_bm[j * N + k] = J_bm_w[(j-3) * N + k];
			}
        
    }

	J_Base2EE(r_e, r_b, J_bE);
	MatrixMulti_(6, 6, N, J_bE, J_bm, J_bE_multi_J_bm);
	MatrixAdd_( 6, N, J_bE_multi_J_bm, Jm, J_g);		
			
	for(int i=0;i<3;i++)
    {
		for(int j=0;j<N;j++)
		{
			J_g_v[i*N+j] = J_g[i*N+j];
			J_g_w[i*N+j] = J_g[(i+3)*N+j];
		}
    }
// ***********************************************************************************************************************************************************************
    
   
    for(double tau = 0; tau < 1; tau += delta_tau)
    {
        for(int i = 0; i < N; i++)
        {
            q_dot_temp = calc_q_dot(tau);
            
            q_dot[i] = q_dot_temp * (para[i] - q_INITIAL[i]);
            
            q[i] = q[i] + q_dot[i] * delta_tau;
        }
        
// ************************** eta_b_dot = (- xi_b.T @ J_bm_w @ q_dot) / 2    ********************************************************************************************
        for(int i = 0; i < 3; i++)
        {
            xi_b_tempt[i] = (-1)*xi_b[i];
        }   
        MatrixMulti_( 1, 3, N, xi_b_tempt, J_bm_w, eta_b_dot_tempt_1);
        MatrixMulti_( 1, N, 1, eta_b_dot_tempt_1, q_dot, eta_b_dot_tempt_2);
        eta_b_dot = (*eta_b_dot_tempt_2) / 2;
//************************************************************************************************************************************************************************
        

// ****************************************** xi_b_dot = ((eta_b * np.eye(3) - cross(xi_b)) @ J_bm_w @ q_dot) / 2  ********************************************************
      
        for(int i = 0; i < 9; i++)
        {
            xi_b_dot_tempt_1[i] = (*eta_b) * eye[i];
        }
        cross(xi_b, cross_xi_b);
        for(int i = 0; i < 9; i++)
        {
            xi_b_dot_tempt_1[i] = xi_b_dot_tempt_1[i] - cross_xi_b[i];
        }
        MatrixMulti_( 3, 3, N, xi_b_dot_tempt_1, J_bm_w, xi_b_dot_tempt_2);
        MatrixMulti_( 3, N, 1, xi_b_dot_tempt_2, q_dot, xi_b_dot_tempt_3);     
        for(int i=0;i<3;i++)
        {
            xi_b_dot[i] = xi_b_dot_tempt_3[i] / 2;
        }
//**********************************************************************************************************************************************************************        

// ************************************** eta_end_dot = (- xi_end.T @ J_g_w @ q_dot) / 2  *******************************************************************************
        for(int i = 0; i < 3; i++)
        {
            xi_end_tempt[i] = (-1)*xi_end[i];
        }
        MatrixMulti_( 1, 3, N, xi_end_tempt, J_g_w, eta_end_dot_tempt_1);
        MatrixMulti_( 1, N, 1, eta_end_dot_tempt_1, q_dot, eta_end_dot_tempt_2);
        eta_end_dot = (*eta_end_dot_tempt_2) / 2;
//**********************************************************************************************************************************************************************        
        
// ************************** xi_end_dot = ((eta_end * np.eye(3) - cross(xi_end)) @ J_g_w @ q_dot) / 2  *****************************************************************
        
        for(int i = 0; i < 9; i++)
        {
            xi_end_dot_tempt_1[i] = (*eta_end) * eye[i];
        }
        cross(xi_end, cross_xi_end);
        for(int i = 0; i < 9; i++)
        {
            xi_end_dot_tempt_1[i] = xi_end_dot_tempt_1[i] - cross_xi_end[i];
        }
        MatrixMulti_( 3, 3, N, xi_end_dot_tempt_1, J_g_w, xi_end_dot_tempt_2);
        MatrixMulti_( 3, N, 1, xi_end_dot_tempt_2, q_dot, xi_end_dot_tempt_3);
        for(int i=0;i<3;i++)
        {
            xi_end_dot[i] = xi_end_dot_tempt_3[i] / 2;
        }
//****************************************************************************************************************************************************************************        
        
// ************************************  next eta_b, xi_b, xi_end *****************************************************************************************************************        
        (*eta_b) = (*eta_b) + eta_b_dot * delta_tau;
        
        for(int i=0;i<3;i++)
            xi_b[i] = xi_b[i] + xi_b_dot[i] * delta_tau;
        
        (*eta_end) = (*eta_end) + eta_end_dot * delta_tau;
        for(int i=0;i<3;i++)
            xi_end[i] = xi_end[i] + xi_end_dot[i] * delta_tau;
        
        MatrixMulti_( 3, N, 1, J_g_v, q_dot, v_e);
        
        for(int i=0;i<3;i++)
            Pe[i] = Pe[i] + v_e[i] * delta_tau;
        
        MatrixMulti_( 3, N, 1, J_bm_v, q_dot, v_b);
        
        for(int i=0;i<3;i++)
            r_b[i] = r_b[i] + v_b[i] * delta_tau;

        quaternion2dc((*eta_b), xi_b[0], xi_b[1], xi_b[2], A_b);      
        links_transform(A_b, q, A_links_transform);       
        calc_r(r_b, A_b, A_links_transform, r);
        
        for(int i=0;i<3;i++)
	    {
		    b_tempt_N_minus_1[i] = b[(N-1)*3+i];		
		    for(int j=0;j<3;j++)
			    A_links_transform_N_minus_1[i*3+j] = A_links_transform[(N-1)*9+i*3+j];
            
	    }
	
	    MatrixMulti_(3, 3, 1, A_links_transform_N_minus_1, b_tempt_N_minus_1, r_e_tempt);

        //  计算末端位置
	    for(int i=0;i<3;i++)
		    r_e[i] = r[(N-1)*3+i] + r_e_tempt[i]; 
        
        
        // 计算 关节位置
	    calc_p(r, A_links_transform, p);
        
        for(int i=0;i<N;i++)
	    {
		        
//*****************************************************************************************************************************************************************
		    MatrixExtract_( 1, 3*N, 1, 1, i*3+1, i*3+3, p, p_tempt);		
		    ScaleMatrix_( 1, 3, (-1), p_tempt, p_tempt);		
		    MatrixAdd_( 1, 3, r_e, p_tempt, r_e_minus_p_i);		
		    MatrixExtract_( 1, 3*3*N, 1, 1, i*9+1, i*9+9, A_links_transform, A_links_transform_tempt_i);
//**************************************************************************************************************************************************************
		    
		    
		    MatrixMulti_(3, 3, 1, A_links_transform_tempt_i, Ez, A_links_transform_tempt_i_multi_Ez);		    
		    cross(A_links_transform_tempt_i_multi_Ez, A_links_transform_tempt_i_multi_Ez_cross);		    
		    MatrixMulti_(3, 3, 1, A_links_transform_tempt_i_multi_Ez_cross, r_e_minus_p_i, Jm_v_tempt_i);
		    
		    for(int k = 0; k < 3; k++)
		    {
			    Jm_v_tempt[i * 3 + k] = Jm_v_tempt_i[k];
			    Jm_w_tempt[i * 3 + k] = A_links_transform_tempt_i_multi_Ez[k];
		    }
		
		    for(int jj = 0; jj < 6; jj++)
		    {
			    if(jj<3)
				    Jm_tempt[i * 6 + jj] = Jm_v_tempt[i * 3 + jj];
			    else
				    Jm_tempt[i * 6 + jj] = Jm_w_tempt[i * 3 + (jj - 3)];
		    }
		
	    }
	
	    MatrixTranspose_(N, 3, Jm_v_tempt, Jm_v);
	    MatrixTranspose_(N, 3, Jm_w_tempt, Jm_w);
	    MatrixTranspose_(N, 6, Jm_tempt, Jm);
	    calc_J_bm(r, r_b, A_b, A_links_transform, p, J_bm_w, J_bm_v);

	    for(int j = 0; j < 6; j++)
        {
		    for(int k = 0; k < N; k++)
			    {
				    if(j < 3)
					    J_bm[j * N + k] = J_bm_v[j * N + k];
				    else
					    J_bm[j * N + k] = J_bm_w[(j-3) * N + k];
			    }
        }

	    J_Base2EE(r_e, r_b, J_bE);
	    MatrixMulti_(6, 6, N, J_bE, J_bm, J_bE_multi_J_bm);
	    for(int i = 0; i < 6; i++)
		    for(int j = 0; j < N; j++)
			    J_g[i * N + j] = J_bE_multi_J_bm[i * N + j] + Jm[i * N + j];
	
	    for(int i=0;i<3;i++)
        {
		    for(int j=0;j<N;j++)
		    {
			    J_g_v[i*N+j] = J_g[i*N+j];
			    J_g_w[i*N+j] = J_g[(i+3)*N+j];
		    }
        } 
        
    }

    
    delete[] q;
    delete[] q_dot;
    delete[] quaternion_base;
    delete[] A_b;
    delete[] xi_b_tempt;

    delete[] xi_end_tempt;
    delete[] A_links_transform;
    delete[] A_b_multi_b_b;
    delete[] r_b_tempt_1;
    delete[] r_b_tempt_2;
    delete[] r_b_tempt_2_temp;
    delete[] a_tempt_i;
    delete[] r_b_tempt_3;
    delete[] r_b;
    delete[] a_tempt_k;
    delete[] b_tempt_k;
    delete[] l_tempt_k;
    delete[] A_links_transform_tempt_i;
    delete[] A_links_transform_tempt_i_multi_a_tempt_i;
    delete[] A_links_transform_tempt_k;
    delete[] A_links_transform_tempt_k_multi_l_tempt_k;
    delete[] r_b_tempt_3_tempt;
    delete[] r;
    delete[] r_e;
    delete[] Pe_initial;
    delete[] A_links_transform_N_minus_1;
    delete[] r_e_tempt;
    delete[] b_tempt_N_minus_1;
    delete[] p;
    delete[] A_links_transform_tempt_i_multi_Ez;
    delete[] A_links_transform_tempt_i_multi_Ez_cross;
    delete[] Jm_v ;
    delete[] Jm_v_tempt;
    delete[] Jm_v_tempt_i ;
    delete[] Jm_w ;
    delete[] Jm_w_tempt;
    delete[] p_tempt;
    delete[] r_e_minus_p_i;
    delete[] Jm ;
    delete[] Jm_tempt;
    delete[] J_bm_w ;
    delete[] J_bm_v;
    delete[] J_bm ;	
    delete[] J_bE ;
    delete[] J_g ;
    delete[] J_bE_multi_J_bm ;
    delete[] J_g_v ;
    delete[] J_g_w ;
    delete[] quaternion_end;
    delete[] eta_b_dot_tempt_1;
    delete[] eta_b_dot_tempt_2;
    delete[] xi_b_dot ;
    delete[] cross_xi_b ;
    delete[] xi_b_dot_tempt_1 ;
    delete[] xi_b_dot_tempt_2 ;
    delete[] xi_b_dot_tempt_3 ;
    delete[] eta_end_dot_tempt_1;
    delete[] eta_end_dot_tempt_2;
    delete[] xi_end_dot ;
    delete[] cross_xi_end ;
    delete[] xi_end_dot_tempt_1 ;
    delete[] xi_end_dot_tempt_2 ;
    delete[] xi_end_dot_tempt_3 ;
    delete[] v_e ;
    delete[] v_b ;
    
    
}
























