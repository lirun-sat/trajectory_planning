
#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void calc_J_bm(double* r, double* r_b, double* A_b, double* A_links_transform, double* p, double* J_bm_w, double* J_bm_v)
{
	/*
    :param I_b_body: 在基座本体系中表示的基座惯量
    :param I_links_body: 在各个连杆本体系中表示的连杆惯量
	*/
	

	// double I_b[9];
    double* I_b;
    I_b = new double[9];
    
    double* I_b_tempt_2;
    double* I_b_tempt_1;
    I_b_tempt_2 = new double[9];
    I_b_tempt_1 = new double[9];
	
    double* r_g;
    r_g = new double[3];
	
	double total_mass = m_b;

	double* a_tempt_1;
	double* a_tempt_2;
	double* A_links_transform_tempt;
	double* A_links_transform_tempt_2;
    // p_tempt_1 = new double[3];
    // p_tempt_2 = new double[3];
    a_tempt_1 = new double[3];
    a_tempt_2 = new double[3];
    A_links_transform_tempt = new double[9];
    A_links_transform_tempt_2 = new double[9];
    
	
	double* I_links_body_tempt;
	double* I_links_tempt_1;
	double* I_links_tempt_2;
	double* I_links;
	double* I_links_i_tempt;
    I_links_body_tempt = new double[9];
    I_links_tempt_1 = new double[9];
    I_links_tempt_2 = new double[9];
    I_links = new double[N * 3 * 3];
    I_links_i_tempt = new double[9];
	
    double* Hw;
    Hw = new double[9];
	
 
    double* r_i_b;
    double* r_i_b_cross_1;
    double* r_i_b_cross_1_T;
    double* r_i_b_cross_1_tempt;
    double* r_i_b_cross_2;
    r_i_b = new double[3];
    r_i_b_cross_1 = new double[9];
    r_i_b_cross_1_T = new double[9];
	r_i_b_cross_1_tempt = new double[9];
    r_i_b_cross_2 = new double[9];
	
    double* joint_revolute_axis_1;
    double* joint_revolute_axis_2;
    joint_revolute_axis_1 = new double[3];
    joint_revolute_axis_2 = new double[3];
	

    double* JR;
    double* JR_i_tempt;
    JR = new double[N * N * 3];
    JR_i_tempt = new double[3 * N];
    
    for(int i = 0; i < N; i++)
        for(int j = 0; j < 3; j++)
            for(int k = 0; k < N; k++)
                JR[i * 3*N + j * N + k] = 0;  // ******************************************************************************************************************************



    double* JT;
    double* JT_i_tempt;
    JT = new double[N * N * 3];
    JT_i_tempt = new double[3 * N];
    
    
    for(int i = 0; i < N; i++)
        for(int j = 0; j < 3; j++)
            for(int k = 0; k < N; k++)
                JT[i * 3*N + j * N + k] = 0;

 
    double* J_Tw;
    double* J_Tw_tempt;
    J_Tw = new double[3*N];
    J_Tw_tempt = new double[3*N];
    
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < N ; j++)
            J_Tw[i * N + j] = 0;
	

    double* J_bm_w_tempt;

    J_bm_w_tempt = new double[3*N];
	
    double* joint_revolute_axis_cross_r_i_minus_p_i;
    double* joint_revolute_axis_cross_r_i_minus_p_j;
    joint_revolute_axis_cross_r_i_minus_p_i = new double[3];
    joint_revolute_axis_cross_r_i_minus_p_j = new double[3];

    double* joint_revolute_axis_2_cross;

    joint_revolute_axis_2_cross = new double[3*3];
	
    double* r_i_minus_p_j;
    double* r_g_minus_r_b;
    double* r_g_minus_r_b_cross;
    double* r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross;
    r_i_minus_p_j = new double[3];
    r_g_minus_r_b = new double[3];
    r_g_minus_r_b_cross = new double[9];
    r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross = new double[9];
    
    double* r_g_minus_r_b_cross_multi_J_Tw;
    r_g_minus_r_b_cross_multi_J_Tw = new double[3*N];
	
    double* H_wq;
    double* H_wq_tempt_1;
    double* H_wq_tempt_2;
    double* H_wq_tempt_3;
    H_wq = new double[3*N];
    H_wq_tempt_1 = new double[3*N];
    H_wq_tempt_2 = new double[3*N];
    H_wq_tempt_3 = new double[3*N];
    
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < N ; j++)
            H_wq[i * N + j] = 0;
    
	
    double* H_s;
    H_s = new double[3*3];

	
	double* H_q;
    H_q = new double[3*N];
	
    double* r_g_minus_r_b_cross_multi_H_s_inv;
    double* r_g_minus_r_b_cross_multi_H_s_inv_multi_H_q;
    r_g_minus_r_b_cross_multi_H_s_inv = new double[3*3];
    r_g_minus_r_b_cross_multi_H_s_inv_multi_H_q = new double[3*N];
    
    double* H_s_inv;
    double* H_s_inv_tempt;
    H_s_inv = new double[3*3];
    H_s_inv_tempt = new double[3*3];
	
	
	// 计算惯性系中表示的基座惯量
    
	MatrixMultiTranspose(I_b_body, A_b, I_b_tempt_1);
	
	MatrixMulti_(3, 3, 3, A_b, I_b_tempt_1, I_b_tempt_2);
	
	MatrixCopy_( 3, 3, I_b_tempt_2, I_b );
	
	// 计算系统质心位置
	
	ScaleMatrix_( 1, 3, m_b, r_b, r_g);
	
	
	for(int i=0;i<N;i++)
	{
		total_mass += m[i];  //  total_mass 初始化为 m_b
		
		for(int j=0;j<3;j++)
			r_g[j] = r_g[j] + m[i] * r[i * 3 + j];
	}
	
	ScaleMatrix_( 1, 3, 1/total_mass, r_g, r_g);
	
    MatrixSub_( 1, 3, r_g, r_b, r_g_minus_r_b);
	
	cross(r_g_minus_r_b, r_g_minus_r_b_cross);
	
	
	/*   // OK
	//-----------------------------------------------------------------------------------------------------------------------------------
	cout << "r_g_minus_r_b_cross:" << "  " 
	               << r_g_minus_r_b_cross[0] << "  "
	               << r_g_minus_r_b_cross[1] << "  "
	               << r_g_minus_r_b_cross[2] << endl;
	               cout << r_g_minus_r_b_cross[3] << "  "
	               << r_g_minus_r_b_cross[4] << "  "
	               << r_g_minus_r_b_cross[5] << endl;
	               cout << r_g_minus_r_b_cross[6] << "  "
	               << r_g_minus_r_b_cross[7] << "  "
	               << r_g_minus_r_b_cross[8] << endl;
	               
	//------------------------------------------------------------------------------------------------------------------------------------
	*/
	
	

	
	
	//  初始化 Hw
	for(int j=0;j<3;j++)
		for(int k=0;k<3;k++)
			Hw[j * 3 + k] = I_b[j * 3 + k];
	
	
	for(int i = 0; i < N; i++)
	{
		for(int j=0;j<3;j++)
			r_i_b[j] = r[i * 3 + j] - r_b[j];
			
		cross(r_i_b, r_i_b_cross_1);  // r_i_b cross
		
		for(int j=0;j<3;j++)
        {
			for(int k=0;k<3;k++)
			{
				r_i_b_cross_1_tempt[j * 3 + k] = r_i_b_cross_1[j * 3 + k];
				
				A_links_transform_tempt[j * 3 + k] = A_links_transform[i * 9 + j * 3 + k];  // A_links_transform_tempt  第 i 个连杆的旋转矩阵
				
				I_links_body_tempt[j * 3 + k] = I_links_body[i * 9 + j * 3 + k];
			}
        }
		
		MatrixTranspose_(3, 3, r_i_b_cross_1, r_i_b_cross_1_T);    // MatrixTranspose(r_i_b_cross_1);  // (r_i -r_b)叉 转置
		MatrixMulti_(3, 3, 3, r_i_b_cross_1_T, r_i_b_cross_1_tempt, r_i_b_cross_2);  //  r_i_b_cross_2 = (r_i -r_b)叉.T * (r_i -r_b)叉
	
		MatrixMulti_(3, 3, 3, A_links_transform_tempt, I_links_body_tempt, I_links_tempt_1);
		
		MatrixMultiTranspose(I_links_tempt_1, A_links_transform_tempt, I_links_tempt_2);  // I_links_tempt_2   第 i 个连杆 在 惯性系 下的 惯量
		
		
		for(int j=0;j<3;j++)
		{	
			for(int k=0;k<3;k++)
		    {
				I_links[i * 9 + j * 3 + k] = I_links_tempt_2[j * 3 + k];
				
				I_links_i_tempt[j * 3 + k] = I_links[i * 9 + j * 3 + k];  //  缓存 第 i 个连杆的 惯量
		
	            /*
                    Hw = I_b
                    for i in range(N):
                        Hw += I_links[i, :, :] + m[i] * cross(r[i, :, :] - r_b).T @ cross(r[i, :, :] - r_b)
                */
				
				Hw[j * 3 + k] += I_links_i_tempt[j * 3 + k] + m[i] * r_i_b_cross_2[j * 3 + k];
		
			}
		}
		
		MatrixMulti_(3, 3, 1, A_links_transform_tempt, Ez, joint_revolute_axis_1);  //  joint_revolute_axis_1  第 i 个关节的旋转轴
		
		for(int j = 0; j < (i+1); j++)
		{
			for(int jj = 0; jj < 3; jj++)
			{
				r_i_minus_p_j[jj] = r[i * 3 + jj] - p[j * 3 + jj];    // 求 r_i - p_j （r_i_minus_p_j） 
				
				for(int kk = 0; kk < 3; kk++)
					A_links_transform_tempt_2[jj * 3 + kk] = A_links_transform[j * 9 + jj * 3 + kk];  //  缓存 第 j 个连杆的旋转矩阵
				
			}
			
			MatrixMulti_(3, 3, 1, A_links_transform_tempt_2, Ez, joint_revolute_axis_2);  // joint_revolute_axis_2  第 j 个关节的旋转轴
			
            cross(joint_revolute_axis_2, joint_revolute_axis_2_cross);				

			MatrixMulti_(3, 3, 1, joint_revolute_axis_2_cross, r_i_minus_p_j, joint_revolute_axis_cross_r_i_minus_p_j);
			
			for(int k = 0; k < 3; k++)
			{	
				if(j == i)
				{
					JR[i * (3*N) + k * N + j] = joint_revolute_axis_1[k];
					JT[i * (3*N) + k * N + j] = joint_revolute_axis_cross_r_i_minus_p_j[k];
				}
				else
				{
					JR[i * (3*N) + k * N + j] = JR[(i-1) * (3*N) + k * N + j];					
					JT[i * (3*N) + k * N + j] = joint_revolute_axis_cross_r_i_minus_p_j[k];			
				}
			}
            
		}
		
		for(int j = 0; j < 3; j++)
        {
			for(int k = 0; k < N; k++)
			{
				JR_i_tempt[j * N + k] = JR[i * (3 * N) + j * N + k];  //  缓存 JR_i
				JT_i_tempt[j * N + k] = JT[i * (3 * N) + j * N + k];  //  缓存 JT_i
				
			}
        }
		
		MatrixMulti_(3, 3, N, I_links_i_tempt, JR_i_tempt, H_wq_tempt_1);
		
		MatrixMulti_(3, 3, N, r_i_b_cross_1_tempt, JT_i_tempt, H_wq_tempt_2);
				
		ScaleMatrix_( 3, N, m[i], H_wq_tempt_2, H_wq_tempt_2);		
		
		MatrixAdd_( 3, N, H_wq_tempt_1, H_wq_tempt_2, H_wq_tempt_3);	
			
		MatrixAdd_( 3, N, H_wq_tempt_3, H_wq, H_wq);
		
		ScaleMatrix_( 3, N, m[i], JT_i_tempt, JT_i_tempt);
		
		MatrixAdd_( 3, N, JT_i_tempt, J_Tw, J_Tw);
		
		
	}
	   
	/*    // OK
	//-----------------------------------------------------------------------------------------------------------------------------------
	cout << "Hw:" << "  " 
	              << Hw[0] << "  "
	              << Hw[1] << "  "
	              << Hw[2] << endl;
	         cout << Hw[3] << "  "
	              << Hw[4] << "  "
	              << Hw[5] << endl;
	         cout << Hw[6] << "  "
	              << Hw[7] << "  "
	              << Hw[8] << endl;
	               
	//------------------------------------------------------------------------------------------------------------------------------------
	*/


	
	
	
	
	MatrixMulti_(3, 3, 3, r_g_minus_r_b_cross, r_g_minus_r_b_cross, r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross);
			
	ScaleMatrix_( 3, 3, total_mass, r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross, r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross);
	
	MatrixAdd_( 3, 3, r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross, Hw, H_s);		
	
	
	/*
	//-----------------------------------------------------------------------------------------------------------------------------------
	cout << "H_s:" << "  " 
	              << H_s[0] << "  "
	              << H_s[1] << "  "
	              << H_s[2] << endl;
	         cout << H_s[3] << "  "
	              << H_s[4] << "  "
	              << H_s[5] << endl;
	         cout << H_s[6] << "  "
	              << H_s[7] << "  "
	              << H_s[8] << endl;
	               
	//------------------------------------------------------------------------------------------------------------------------------------
	*/
	
	
	
	
	// Hq = Hwq - cross(r_g - r_b) @ JTw
	MatrixMulti_(3, 3, N, r_g_minus_r_b_cross, J_Tw, r_g_minus_r_b_cross_multi_J_Tw);
	
	/*
	for(int kk = 0; kk < 3; kk++)
		for(int jj = 0; jj < N; jj++)
			H_q[kk * N + jj] = H_wq[kk * N + jj] - r_g_minus_r_b_cross_multi_J_Tw[kk * N + jj];
	
	*/		
			
	MatrixSub_( 3, N, H_wq, r_g_minus_r_b_cross_multi_J_Tw, H_q);
    
    
    /*    // OK
	//-----------------------------------------------------------------------------------------------------------------------------------
	cout << "H_q:" << "  " 
	              << H_q[0] << "  "
	              << H_q[1] << "  "
	              << H_q[2] << "  "
	              << H_q[3] << "  "
	              << H_q[4] << "  "
	              << H_q[5] << "  "
	              << H_q[6] << endl;
	         cout << H_q[7] << "  "
	              << H_q[8] << "  "
	              << H_q[9] << "  "
	              << H_q[10] << "  "
	              << H_q[11] << "  "
	              << H_q[12] << "  "
	              << H_q[13] << endl;
	         cout << H_q[14] << "  "
	              << H_q[15] << "  "
	              << H_q[16] << "  "
	              << H_q[17] << "  "
	              << H_q[18] << "  "
	              << H_q[19] << "  "
	              << H_q[20] << endl;
	               
	//------------------------------------------------------------------------------------------------------------------------------------
	*/
	
	
			    
	//J_bm_w = (-1) * linalg.inv(Hs) @ Hq
    //J_bm_v = (-1) * (JTw / (m_b + sum(m)) + cross(r_g - r_b) @ linalg.inv(Hs) @ Hq)
    
	
	LUP_solve_inverse(H_s, 3, H_s_inv);	
	
	MatrixCopy_( 3, 3, H_s_inv, H_s_inv_tempt );
	
	/*
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			H_s_inv_tempt[i * 3 + j] = H_s_inv[i * 3 + j];
	*/
	
	MatrixMulti_(3, 3, N, H_s_inv_tempt, H_q, J_bm_w_tempt);
	
	ScaleMatrix_( 3, N, (-1), J_bm_w_tempt, J_bm_w);
	
	ScaleMatrix_( 3, N, 1/total_mass, J_Tw, J_Tw_tempt);
	
	MatrixMulti_(3, 3, 3, r_g_minus_r_b_cross, H_s_inv_tempt, r_g_minus_r_b_cross_multi_H_s_inv);
	MatrixMulti_(3, 3, N, r_g_minus_r_b_cross_multi_H_s_inv, H_q, r_g_minus_r_b_cross_multi_H_s_inv_multi_H_q);
	
	MatrixAdd_( 3, N, J_Tw_tempt, r_g_minus_r_b_cross_multi_H_s_inv_multi_H_q, J_bm_v);
	
	ScaleMatrix_( 3, N, (-1), J_bm_v, J_bm_v);
	
		
	/*	
	//-------------------OK--------------------------------------------------------------------------------------------------
	for(int i=0;i<3;i++)
		for(int j=0;j<N;j++)
			cout << J_bm_v[i * N + j] << endl; 
	//-----------------------------------------------------------------------------------------------------------------------
	
	*/
	
	
    
    delete[] I_b;
    delete[] I_b_tempt_2;
    delete[] I_b_tempt_1;
    delete[] r_g;
	delete[] a_tempt_1;
	delete[] a_tempt_2;
	delete[] A_links_transform_tempt;
	delete[] A_links_transform_tempt_2;
	delete[] I_links_body_tempt;
	delete[] I_links_tempt_1;
	delete[] I_links_tempt_2;
	delete[] I_links;
	delete[] I_links_i_tempt;
    delete[] Hw;	
    delete[] r_i_b;
    delete[] r_i_b_cross_1;
    delete[] r_i_b_cross_1_T;
    delete[] r_i_b_cross_1_tempt;
    delete[] r_i_b_cross_2;
    delete[] joint_revolute_axis_1;
    delete[] joint_revolute_axis_2;
    delete[] JR;
    delete[] JR_i_tempt;
    delete[] JT;
    delete[] JT_i_tempt;
    delete[] J_Tw;
    delete[] J_Tw_tempt;
    delete[] J_bm_w_tempt;
    delete[] joint_revolute_axis_cross_r_i_minus_p_i;
    delete[] joint_revolute_axis_cross_r_i_minus_p_j;
    delete[] joint_revolute_axis_2_cross;
    delete[] r_i_minus_p_j;
    delete[] r_g_minus_r_b;
    delete[] r_g_minus_r_b_cross;
    delete[] r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross;
    delete[] r_g_minus_r_b_cross_multi_J_Tw;
    delete[] H_wq;
    delete[] H_wq_tempt_1;
    delete[] H_wq_tempt_2;
    delete[] H_wq_tempt_3;
    delete[] H_s;
	delete[] H_q;	
    delete[] r_g_minus_r_b_cross_multi_H_s_inv;
    delete[] r_g_minus_r_b_cross_multi_H_s_inv_multi_H_q;
    delete[] H_s_inv;
    delete[] H_s_inv_tempt;
	
}


























