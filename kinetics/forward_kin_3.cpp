#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include "data.h"
#include "utils.h"
#include "openGJK.h"


void forward_kin_3(double* para, double* eta_end, double* xi_end, double* Pe, double* eta_b, double* xi_b, 
                   double* p_e_initial, double* locus, double* delta_xi_b_distrb_max, double* manipl, double* T_min, double* collision)
{
    double* q                                         = new double[N];
    double* q_dot                                     = new double[N];
    double* q_ddot                                    = new double[N];
    double* quaternion_base_initial                   = new double[4];
    double* A_b                                       = new double[3*3];
    double* xi_b_tempt                                = new double[3];
    double* xi_end_tempt                              = new double[3];
    double* A_links_transform                         = new double[3*3*N];
    double* A_b_multi_b_b                             = new double[3];
    double* r_b_tempt_1                               = new double[3];
    double* r_b_tempt_2                               = new double[3];
    double* r_b_tempt_2_temp                          = new double[3];
    double* a_tempt_i                                 = new double[3];
    double* r_b_tempt_3                               = new double[3];
    double* r_b                                       = new double[3];
    double* a_tempt_k                                 = new double[3];
    double* b_tempt_k                                 = new double[3];
    double* l_tempt_k                                 = new double[3];    
    double* A_links_transform_tempt_i                 = new double[3*3];
    double* A_links_transform_tempt_i_multi_a_tempt_i = new double[3];
    double* A_links_transform_tempt_k                 = new double[3*3];
    double* A_links_transform_tempt_k_multi_l_tempt_k = new double[3];
    double* r_b_tempt_3_tempt                         = new double[3];
    double* r                                         = new double[3 * N];
    double* r_e                                       = new double[3];
    double* Pe_initial                                = new double[3];
    double* A_links_transform_N_minus_1               = new double[3*3];
    double* r_e_tempt                                 = new double[3];
    double* b_tempt_N_minus_1                         = new double[3];
    double* p                                         = new double[N*3];
    double* A_links_transform_tempt_i_multi_Ez        = new double[3];
    double* A_links_transform_tempt_i_multi_Ez_cross  = new double[9];
    double* Jm_v                                      = new double[3*N];
    double* Jm_v_tempt                                = new double[N*3];
    double* Jm_v_tempt_i                              = new double[3];
    double* Jm_w                                      = new double[3*N];
    double* Jm_w_tempt                                = new double[N*3];
    double* p_tempt                                   = new double[3];
    double* r_e_minus_p_i                             = new double[3];
    double* Jm                                        = new double[6*N];
    double* Jm_tempt                                  = new double[N*6];
    double* J_bm_w                                    = new double[3*N];
    double* J_bm_v                                    = new double[3*N];
    double* J_bm                                      = new double[6*N];	
    double* J_bE                                      = new double[6*6];
    double* J_g                                       = new double[6*N];
    double* J_bE_multi_J_bm                           = new double[6*N];
    double* J_g_v                                     = new double[3*N];
    double* J_g_w                                     = new double[3*N];    
    double* quaternion_end                            = new double[4];
    double* eta_b_dot_tempt_1                         = new double[N];
    double* eta_b_dot_tempt_2                         = new double;  
    double* xi_b_dot                                  = new double[3];
    double* cross_xi_b                                = new double[3*3];  
    double* xi_b_dot_tempt_1                          = new double[3*3];
    double* xi_b_dot_tempt_2                          = new double[3*N];
    double* xi_b_dot_tempt_3                          = new double[3]; 
    double* eta_end_dot_tempt_1                       = new double[N];
    double* eta_end_dot_tempt_2                       = new double;
    double* xi_end_dot                                = new double[3];
    double* cross_xi_end                              = new double[3*3]; 
    double* xi_end_dot_tempt_1                        = new double[3*3];
    double* xi_end_dot_tempt_2                        = new double[3*N];
    double* xi_end_dot_tempt_3                        = new double[3];
    double* v_e                                       = new double[3];
    double* v_b                                       = new double[3];
    double* A_b_expand                                = new double[2*3 * 2*3];
	double* A_b_expand_transpose                      = new double[2*3 * 2*3];
	double* A_b_expand_transpose_multi_J_g            = new double[6 * N];
	double* A_b_expand_transpose_multi_J_g_transpose  = new double[N * 6];
	double* manipl_temp                               = new double[6 * 6];
    double* xi_b_initial                              = new double[3];
    double* delta_xi_base                             = new double[3];
    double* delta_eta_base_tempt                      = new double;
    double* delta_eta_base                            = new double;
    double* delta_xi_base_tempt                       = new double[3];

    double eta_end_dot = 0;
    double eta_b_dot = 0;
    double total_mass_for_links = 0;
	double total_mass = 0;	
	double q_dot_temp = 0;
    double q_ddot_temp = 0;
    double eta_b_initial = 0;
    double delta_xi_base_norm_max = 0;
    double delta_xi_base_norm = 0;
    double q_dot_max = 0;
    double q_ddot_max = 0;
    double q_dot_max_vmax = 0;
    double q_ddot_max_alphamax = 0;
    double T_min_temp = 0;

    double* base_object_tempt          = new double[nvrtx * 3];
    double* left_solar_object_tempt    = new double[nvrtx * 3];
    double* right_solar_object_tempt   = new double[nvrtx * 3];
    double* link_1_part_1_object_tempt = new double[nvrtx * 3];
    double* link_1_part_2_object_tempt = new double[nvrtx * 3];
    double* link_2_part_1_object_tempt = new double[nvrtx * 3];
    double* link_2_part_2_object_tempt = new double[nvrtx * 3];
    double* link_5_part_1_object_tempt = new double[nvrtx * 3];
    double* link_5_part_2_object_tempt = new double[nvrtx * 3];
    double* link_5_part_3_object_tempt = new double[nvrtx * 3];
    double* link_5_part_4_object_tempt = new double[nvrtx * 3];
    double* link_5_part_5_object_tempt = new double[nvrtx * 3];
    double* link_6_part_1_object_tempt = new double[nvrtx * 3];
    double* link_6_part_2_object_tempt = new double[nvrtx * 3];
    double* link_7_object_tempt        = new double[nvrtx * 3];

    double* A_links_transform_0 = new double[3*3];
    double* A_links_transform_1 = new double[3*3];
    double* A_links_transform_4 = new double[3*3];
    double* A_links_transform_5 = new double[3*3];
    double* A_links_transform_6 = new double[3*3];

    double** base_object          = new double*[nvrtx];
    double** left_solar_object    = new double*[nvrtx];
    double** right_solar_object   = new double*[nvrtx];
    double** link_1_part_1_object = new double*[nvrtx];
    double** link_1_part_2_object = new double*[nvrtx];
    double** link_2_part_1_object = new double*[nvrtx];
    double** link_2_part_2_object = new double*[nvrtx];
    double** link_5_part_1_object = new double*[nvrtx];
    double** link_5_part_2_object = new double*[nvrtx];
    double** link_5_part_3_object = new double*[nvrtx];
    double** link_5_part_4_object = new double*[nvrtx];
    double** link_5_part_5_object = new double*[nvrtx];
    double** link_6_part_1_object = new double*[nvrtx];
    double** link_6_part_2_object = new double*[nvrtx];
    double** link_7_object        = new double*[nvrtx];
    
    for(int i = 0; i < nvrtx; i++){
        base_object[i]          = new double[3];
        left_solar_object[i]    = new double[3];
        right_solar_object[i]   = new double[3];
        link_1_part_1_object[i] = new double[3];
        link_1_part_2_object[i] = new double[3];
        link_2_part_1_object[i] = new double[3];
        link_2_part_2_object[i] = new double[3];
        link_5_part_1_object[i] = new double[3];
        link_5_part_2_object[i] = new double[3];
        link_5_part_3_object[i] = new double[3];
        link_5_part_4_object[i] = new double[3];
        link_5_part_5_object[i] = new double[3];
        link_6_part_1_object[i] = new double[3];
        link_6_part_2_object[i] = new double[3];
        link_7_object[i]        = new double[3];
    }

    double* base_center2vertex_temp                 = new double[nvrtx * 3]; 
    double* base_center2left_solar_vertex_temp      = new double[nvrtx * 3]; 
    double* base_center2right_solar_vertex_temp     = new double[nvrtx * 3]; 
    double* link_1_center2link_1_part_1_vertex_temp = new double[nvrtx * 3]; 
    double* link_1_center2link_1_part_2_vertex_temp = new double[nvrtx * 3]; 
    double* link_2_center2link_2_part_1_vertex_temp = new double[nvrtx * 3]; 
    double* link_2_center2link_2_part_2_vertex_temp = new double[nvrtx * 3]; 
    double* link_5_center2link_5_part_1_vertex_temp = new double[nvrtx * 3]; 
    double* link_5_center2link_5_part_2_vertex_temp = new double[nvrtx * 3]; 
    double* link_5_center2link_5_part_3_vertex_temp = new double[nvrtx * 3]; 
    double* link_5_center2link_5_part_4_vertex_temp = new double[nvrtx * 3]; 
    double* link_5_center2link_5_part_5_vertex_temp = new double[nvrtx * 3]; 
    double* link_6_center2link_6_part_1_vertex_temp = new double[nvrtx * 3]; 
    double* link_6_center2link_6_part_2_vertex_temp = new double[nvrtx * 3]; 
    double* link_7_center2link_7_vertex_temp        = new double[nvrtx * 3]; 

    double* distance_base2link_567   = new double[8];
    double* distance_left2link_567   = new double[7];
    double* distance_right2link_567  = new double[7];
    double* distance_link_12link_567 = new double[14];
    double* distance_link_22link_567 = new double[14];


    gkSimplex s;
    s.nvrtx = 0;
    
    gkPolytope bd_base_object;
    gkPolytope bd_left_solar_object;
    gkPolytope bd_right_solar_object;
    gkPolytope bd_link_1_part_1_object;
    gkPolytope bd_link_1_part_2_object;
    gkPolytope bd_link_2_part_1_object;
    gkPolytope bd_link_2_part_2_object;
    gkPolytope bd_link_5_part_1_object;
    gkPolytope bd_link_5_part_2_object;
    gkPolytope bd_link_5_part_3_object;
    gkPolytope bd_link_5_part_4_object;
    gkPolytope bd_link_5_part_5_object;
    gkPolytope bd_link_6_part_1_object;
    gkPolytope bd_link_6_part_2_object;
    gkPolytope bd_link_7_object;
    

    double collision_test_base_link567   = 0;
    double collision_test_left_link567   = 0;
    double collision_test_right_link567  = 0;
    double collision_test_link_1_link567 = 0;
    double collision_test_link_2_link567 = 0;

    double collision_test = 0;

    double distance_base2link_5_part_1          = 0;
    double distance_base2link_5_part_2          = 0;
    double distance_base2link_5_part_3          = 0;
    double distance_base2link_5_part_4          = 0;
    double distance_base2link_5_part_5          = 0;
    double distance_base2link_6_part_1          = 0;
    double distance_base2link_6_part_2          = 0;
    double distance_base2link_7                 = 0;
    double distance_left2link_5_part_2          = 0;
    double distance_left2link_5_part_3          = 0;
    double distance_left2link_5_part_4          = 0;
    double distance_left2link_5_part_5          = 0;
    double distance_left2link_6_part_1          = 0;
    double distance_left2link_6_part_2          = 0;
    double distance_left2link_7                 = 0;
    double distance_right2link_5_part_2         = 0;
    double distance_right2link_5_part_3         = 0;
    double distance_right2link_5_part_4         = 0;
    double distance_right2link_5_part_5         = 0;
    double distance_right2link_6_part_1         = 0;
    double distance_right2link_6_part_2         = 0;
    double distance_right2link_7                = 0;
    double distance_link_1_part_1_link_5_part_2 = 0;
    double distance_link_1_part_1_link_5_part_3 = 0;
    double distance_link_1_part_1_link_5_part_4 = 0;
    double distance_link_1_part_1_link_5_part_5 = 0;
    double distance_link_1_part_1_link_6_part_1 = 0;
    double distance_link_1_part_1_link_6_part_2 = 0;
    double distance_link_1_part_1_link_7        = 0; 
    double distance_link_1_part_2_link_5_part_2 = 0;
    double distance_link_1_part_2_link_5_part_3 = 0;
    double distance_link_1_part_2_link_5_part_4 = 0;
    double distance_link_1_part_2_link_5_part_5 = 0;
    double distance_link_1_part_2_link_6_part_1 = 0;
    double distance_link_1_part_2_link_6_part_2 = 0;
    double distance_link_1_part_2_link_7        = 0;
    double distance_link_2_part_1_link_5_part_2 = 0;
    double distance_link_2_part_1_link_5_part_3 = 0;
    double distance_link_2_part_1_link_5_part_4 = 0;
    double distance_link_2_part_1_link_5_part_5 = 0;
    double distance_link_2_part_1_link_6_part_1 = 0;
    double distance_link_2_part_1_link_6_part_2 = 0;
    double distance_link_2_part_1_link_7        = 0;
    double distance_link_2_part_2_link_5_part_2 = 0;
    double distance_link_2_part_2_link_5_part_3 = 0;
    double distance_link_2_part_2_link_5_part_4 = 0;
    double distance_link_2_part_2_link_5_part_5 = 0;
    double distance_link_2_part_2_link_6_part_1 = 0;
    double distance_link_2_part_2_link_6_part_2 = 0;
    double distance_link_2_part_2_link_7        = 0;

    double distance_limit = 0.01;

    // transform two dimensional array to one array
    
    for(int i = 0; i < nvrtx; i++){
        for(int j = 0; j < 3; j++)
        {
            base_center2vertex_temp[i*3 + j]                 = base_center2vertex[i][j];
            base_center2left_solar_vertex_temp[i*3 + j]      = base_center2left_solar_vertex[i][j];
            base_center2right_solar_vertex_temp[i*3 + j]     = base_center2right_solar_vertex[i][j];
            link_1_center2link_1_part_1_vertex_temp[i*3 + j] = link_1_center2link_1_part_1_vertex[i][j];
            link_1_center2link_1_part_2_vertex_temp[i*3 + j] = link_1_center2link_1_part_2_vertex[i][j];
            link_2_center2link_2_part_1_vertex_temp[i*3 + j] = link_2_center2link_2_part_1_vertex[i][j];
            link_2_center2link_2_part_2_vertex_temp[i*3 + j] = link_2_center2link_2_part_2_vertex[i][j];
            link_5_center2link_5_part_1_vertex_temp[i*3 + j] = link_5_center2link_5_part_1_vertex[i][j];
            link_5_center2link_5_part_2_vertex_temp[i*3 + j] = link_5_center2link_5_part_2_vertex[i][j];
            link_5_center2link_5_part_3_vertex_temp[i*3 + j] = link_5_center2link_5_part_3_vertex[i][j];
            link_5_center2link_5_part_4_vertex_temp[i*3 + j] = link_5_center2link_5_part_4_vertex[i][j];
            link_5_center2link_5_part_5_vertex_temp[i*3 + j] = link_5_center2link_5_part_5_vertex[i][j];
            link_6_center2link_6_part_1_vertex_temp[i*3 + j] = link_6_center2link_6_part_1_vertex[i][j];
            link_6_center2link_6_part_2_vertex_temp[i*3 + j] = link_6_center2link_6_part_2_vertex[i][j];
            link_7_center2link_7_vertex_temp[i*3 + j]        = link_7_center2link_7_vertex[i][j];            
        }
    }

    MatrixTranspose_(nvrtx, 3, base_center2vertex_temp, base_center2vertex_temp);
    MatrixTranspose_(nvrtx, 3, base_center2left_solar_vertex_temp, base_center2left_solar_vertex_temp);
    MatrixTranspose_(nvrtx, 3, base_center2right_solar_vertex_temp, base_center2right_solar_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_1_center2link_1_part_1_vertex_temp, link_1_center2link_1_part_1_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_1_center2link_1_part_2_vertex_temp, link_1_center2link_1_part_2_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_2_center2link_2_part_1_vertex_temp, link_2_center2link_2_part_1_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_2_center2link_2_part_2_vertex_temp, link_2_center2link_2_part_2_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_5_center2link_5_part_1_vertex_temp, link_5_center2link_5_part_1_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_5_center2link_5_part_2_vertex_temp, link_5_center2link_5_part_2_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_5_center2link_5_part_3_vertex_temp, link_5_center2link_5_part_3_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_5_center2link_5_part_4_vertex_temp, link_5_center2link_5_part_4_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_5_center2link_5_part_5_vertex_temp, link_5_center2link_5_part_5_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_6_center2link_6_part_1_vertex_temp, link_6_center2link_6_part_1_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_6_center2link_6_part_2_vertex_temp, link_6_center2link_6_part_2_vertex_temp);
    MatrixTranspose_(nvrtx, 3, link_7_center2link_7_vertex_temp, link_7_center2link_7_vertex_temp);

//**************************************************************************************************************************************************
    for(int i = 0; i < N; i++){
        q[i] = q_INITIAL[i];
    }
        
    zyx2quaternion(RPY_BASE_INITIAL[2], RPY_BASE_INITIAL[1], RPY_BASE_INITIAL[0], quaternion_base_initial);
    quaternion2dc(quaternion_base_initial[0], quaternion_base_initial[1], quaternion_base_initial[2], quaternion_base_initial[3], A_b);

    *eta_b = quaternion_base_initial[0];
    eta_b_initial = quaternion_base_initial[0];

    for(int i=0;i<3;i++){
        xi_b[i] = quaternion_base_initial[i+1];
        xi_b_initial[i] = quaternion_base_initial[i+1];
    }

	*locus = 0;

// ***************************************************************************************************************************************
    links_transform(A_b, q, A_links_transform);
      
    for(int i=0;i<3;i++){
        r_b_tempt_2[i] = 0;
    }
    
    for(int i=0;i<3;i++){
        r_b_tempt_3[i] = 0;
    }
    
    // 自由漂浮，初始动量为 0 ， 以 系统质心为 惯性系原点。
    for (int i = 0; i < N; i++){
        total_mass_for_links += m[i];
    }

	total_mass = m_b + total_mass_for_links;
	MatrixMulti_(3, 3, 1, A_b, b_b, A_b_multi_b_b);	
	ScaleMatrix_( 3, 1, total_mass_for_links, A_b_multi_b_b, r_b_tempt_1);
	
	for(int i=0;i<N;i++){
		for(int j=0;j<3;j++){
			a_tempt_i[j] = a[i*3+j];			
			r_b_tempt_3_tempt[j] = 0;			
			for(int k=0;k<3;k++){
                A_links_transform_tempt_i[j*3+k] = A_links_transform[i*9+j*3+k];
            }
		}
		MatrixMulti_( 3, 3, 1, A_links_transform_tempt_i, a_tempt_i, A_links_transform_tempt_i_multi_a_tempt_i);		
		ScaleMatrix_( 1, 3, m[i], A_links_transform_tempt_i_multi_a_tempt_i, r_b_tempt_2_temp);		
		MatrixAdd_( 1, 3, r_b_tempt_2, r_b_tempt_2_temp, r_b_tempt_2);

		for(int k=0;k<i;k++){
			for(int jj=0;jj<3;jj++){
				a_tempt_k[jj] = a[k*3+jj];
				b_tempt_k[jj] = b[k*3+jj];
				l_tempt_k[jj] = a_tempt_k[jj] + b_tempt_k[jj];				
				for(int kk=0;kk<3;kk++){
                    A_links_transform_tempt_k[jj*3+kk] = A_links_transform[k*9+jj*3+kk];
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
	
    //  计算末端位置
	for(int i=0;i<3;i++){
		b_tempt_N_minus_1[i] = b[(N-1)*3+i];		
		for(int j=0;j<3;j++){
            A_links_transform_N_minus_1[i*3+j] = A_links_transform[(N-1)*9+i*3+j];
        }
			
	}
	MatrixMulti_(3, 3, 1, A_links_transform_N_minus_1, b_tempt_N_minus_1, r_e_tempt);

	for(int i=0;i<3;i++){
        r_e[i] = r[(N-1)*3+i] + r_e_tempt[i];
    }
		
// ********************************************** 初始末端位姿 ******************************************************************************************   	
    for(int i=0;i<3;i++){
        Pe_initial[i] = r_e[i];  
    }
    
    zyx2quaternion(RPY_END_INITIAL[2], RPY_END_INITIAL[1], RPY_END_INITIAL[0], quaternion_end);  
    *eta_end = quaternion_end[0];

    for(int i=0;i<3;i++){
        xi_end[i] = quaternion_end[i+1];
    }
    
    for(int i=0;i<3;i++){
        Pe[i] = Pe_initial[i];
        p_e_initial[i] = Pe_initial[i];
    }
        
// ***********************************************************************************************************************************************************      
    // 计算 关节位置
	calc_p(r, A_links_transform, p);    

// *************************************  计算 初始 Jm_v, Jm_w, Jm, J_bm_w, J_bm_v, J_bm, J_g_v, J_g_w, J_g  ****************************************************
    for(int i=0;i<N;i++){
		MatrixExtract_( 1, 3*N, 1, 1, i*3+1, i*3+3, p, p_tempt);		
		ScaleMatrix_( 1, 3, (-1), p_tempt, p_tempt);		
		MatrixAdd_( 1, 3, r_e, p_tempt, r_e_minus_p_i);		
		MatrixExtract_( 1, 3*3*N, 1, 1, i*9+1, i*9+9, A_links_transform, A_links_transform_tempt_i);						
		MatrixMulti_(3, 3, 1, A_links_transform_tempt_i, Ez, A_links_transform_tempt_i_multi_Ez);		
		cross(A_links_transform_tempt_i_multi_Ez, A_links_transform_tempt_i_multi_Ez_cross);		
		MatrixMulti_(3, 3, 1, A_links_transform_tempt_i_multi_Ez_cross, r_e_minus_p_i, Jm_v_tempt_i);		
		
		for(int k = 0; k < 3; k++){
			Jm_v_tempt[i * 3 + k] = Jm_v_tempt_i[k];
			Jm_w_tempt[i * 3 + k] = A_links_transform_tempt_i_multi_Ez[k];
		}
		
		for(int jj = 0; jj < 6; jj++){
			if(jj<3){
				Jm_tempt[i * 6 + jj] = Jm_v_tempt[i * 3 + jj];
            }
			else{
				Jm_tempt[i * 6 + jj] = Jm_w_tempt[i * 3 + (jj - 3)];
            }
		}
	}
	MatrixTranspose_(N, 3, Jm_v_tempt, Jm_v);
	MatrixTranspose_(N, 3, Jm_w_tempt, Jm_w);
	MatrixTranspose_(N, 6, Jm_tempt, Jm);

    calc_J_bm( r, r_b, A_b, A_links_transform, p, J_bm_w, J_bm_v);

	for(int j = 0; j < 6; j++){
		for(int k = 0; k < N; k++){
			if(j < 3){
				J_bm[j * N + k] = J_bm_v[j * N + k];
            }
			else{
				J_bm[j * N + k] = J_bm_w[(j-3) * N + k];
            }
		}
    }

	J_Base2EE(r_e, r_b, J_bE);
	MatrixMulti_(6, 6, N, J_bE, J_bm, J_bE_multi_J_bm);
	MatrixAdd_( 6, N, J_bE_multi_J_bm, Jm, J_g);		
			
	for(int i=0;i<3;i++){
		for(int j=0;j<N;j++){
			J_g_v[i*N+j] = J_g[i*N+j];
			J_g_w[i*N+j] = J_g[(i+3)*N+j];
		}
    }
// *******************************************  进入迭代计算, 步长 delta_tau  **************************************************************************************
    
    for(double tau = 0; tau < 1; tau += delta_tau){
        for(int i = 0; i < N; i++){
            q_dot_temp = calc_q_dot(tau);
            q_ddot_temp = calc_q_ddot(tau);
            q_dot[i] = q_dot_temp * (para[i] - q_INITIAL[i]);
            q_ddot[i] = q_ddot_temp * (para[i] - q_INITIAL[i]);
            q[i] = q[i] + q_dot[i] * delta_tau;

            if(fabs(q_dot[i]) > q_dot_max){
                q_dot_max = fabs(q_dot[i]);
            }
            
            if(fabs(q_ddot[i]) > q_ddot_max){
                q_ddot_max = fabs(q_ddot[i]);
            }
        }
        if(q_dot_max / joint_angle_velocity_max_limit > q_dot_max_vmax){
            q_dot_max_vmax = q_dot_max / joint_angle_velocity_max_limit;
        }

        if(q_ddot_max / joint_angle_acceleration_max_limit > q_ddot_max_alphamax){
            q_ddot_max_alphamax = q_ddot_max / joint_angle_acceleration_max_limit;
        }

        if(sqrt(q_ddot_max_alphamax > T_min_temp)){
            T_min_temp = sqrt(q_ddot_max_alphamax);
        }
        
        if(q_dot_max_vmax > T_min_temp){
            T_min_temp = q_dot_max_vmax;
        }
        
        // ************************** eta_b_dot = (- xi_b.T @ J_bm_w @ q_dot) / 2    **************************************************************************
        for(int i = 0; i < 3; i++){
            xi_b_tempt[i] = (-1)*xi_b[i];
        }   
        MatrixMulti_( 1, 3, N, xi_b_tempt, J_bm_w, eta_b_dot_tempt_1);
        MatrixMulti_( 1, N, 1, eta_b_dot_tempt_1, q_dot, eta_b_dot_tempt_2);
        eta_b_dot = (*eta_b_dot_tempt_2) / 2;
        //********************************************************************************************************************************************************
        
		// ****************************************** xi_b_dot = ((eta_b * np.eye(3) - cross(xi_b)) @ J_bm_w @ q_dot) / 2  ***************************************
      
        for(int i = 0; i < 9; i++){
            xi_b_dot_tempt_1[i] = (*eta_b) * eye[i];
        }
        cross(xi_b, cross_xi_b);
        for(int i = 0; i < 9; i++){
            xi_b_dot_tempt_1[i] = xi_b_dot_tempt_1[i] - cross_xi_b[i];
        }
        MatrixMulti_( 3, 3, N, xi_b_dot_tempt_1, J_bm_w, xi_b_dot_tempt_2);
        MatrixMulti_( 3, N, 1, xi_b_dot_tempt_2, q_dot, xi_b_dot_tempt_3);     
        for(int i=0;i<3;i++){
            xi_b_dot[i] = xi_b_dot_tempt_3[i] / 2;
        }
		//****************************************************************************************************************************************************        

		// ************************************** eta_end_dot = (- xi_end.T @ J_g_w @ q_dot) / 2  **************************************************************
        for(int i = 0; i < 3; i++){
            xi_end_tempt[i] = (-1)*xi_end[i];
        }
        MatrixMulti_( 1, 3, N, xi_end_tempt, J_g_w, eta_end_dot_tempt_1);
        MatrixMulti_( 1, N, 1, eta_end_dot_tempt_1, q_dot, eta_end_dot_tempt_2);
        eta_end_dot = (*eta_end_dot_tempt_2) / 2;
		//********************************************************************************************************************************************************        
        
		// ************************** xi_end_dot = ((eta_end * np.eye(3) - cross(xi_end)) @ J_g_w @ q_dot) / 2  **************************************************
        
        for(int i = 0; i < 9; i++){
            xi_end_dot_tempt_1[i] = (*eta_end) * eye[i];
        }
        cross(xi_end, cross_xi_end);
        for(int i = 0; i < 9; i++){
            xi_end_dot_tempt_1[i] = xi_end_dot_tempt_1[i] - cross_xi_end[i];
        }
        MatrixMulti_( 3, 3, N, xi_end_dot_tempt_1, J_g_w, xi_end_dot_tempt_2);
        MatrixMulti_( 3, N, 1, xi_end_dot_tempt_2, q_dot, xi_end_dot_tempt_3);
        for(int i=0;i<3;i++){
            xi_end_dot[i] = xi_end_dot_tempt_3[i] / 2;
        }
		//******************************************************************************************************************************************
        
		// ****************  next eta_b, xi_b, eta_end, xi_end, Pe, r_b, A_b, A_links_transform, r, r_e, p, Jm, J_bm, J_g ****************************
        (*eta_b) = (*eta_b) + eta_b_dot * delta_tau;
        for(int i = 0; i < 3; i++){
            xi_b[i] = xi_b[i] + xi_b_dot[i] * delta_tau;
        }
        //**********************************************************  Normalization **************************************************************
        (*eta_b) = (*eta_b) / sqrt((*eta_b) * (*eta_b) + xi_b[0] * xi_b[0] + xi_b[1] * xi_b[1] + xi_b[2] * xi_b[2]);
        for(int i = 0; i < 3; i++){
            xi_b[i] = xi_b[i] / sqrt((*eta_b) * (*eta_b) + xi_b[0] * xi_b[0] + xi_b[1] * xi_b[1] + xi_b[2] * xi_b[2]);
        }
		// **************************************************  计算基座扰动 RPY and delta_xi_b 最大值  *************************************************
    
        cross(xi_b, cross_xi_b);
        MatrixMulti_(3, 3, 1, cross_xi_b, xi_b_initial, delta_xi_base_tempt);
        for(int i = 0; i < 3; i++){
            delta_xi_base_tempt[i] = (-1) * delta_xi_base_tempt[i] - eta_b_initial * xi_b[i];     
            delta_xi_base_tempt[i] = delta_xi_base_tempt[i] + (*eta_b) * xi_b_initial[i];
        }
        for(int i=0;i<3;i++){
            delta_xi_base[i] = delta_xi_base_tempt[i];
        }
    
        MatrixMulti_(1, 3, 1, xi_b, xi_b_initial, delta_eta_base_tempt);
        (*delta_eta_base_tempt) = (*delta_eta_base_tempt) + (*eta_b) * eta_b_initial;
        (*delta_eta_base) = (*delta_eta_base_tempt);

        delta_xi_base_norm = sqrt(delta_xi_base[0] * delta_xi_base[0] + delta_xi_base[1] * delta_xi_base[1] + delta_xi_base[2] * delta_xi_base[2]);

        if(delta_xi_base_norm > delta_xi_base_norm_max){
            delta_xi_base_norm_max = delta_xi_base_norm;
        }

		// ********************************************************************************************************************************
        
        (*eta_end) = (*eta_end) + eta_end_dot * delta_tau;
        for(int i=0;i<3;i++){
            xi_end[i] = xi_end[i] + xi_end_dot[i] * delta_tau;
        }
        //**********************************************************  Normalization ***************************************************
        (*eta_end) = (*eta_end) / sqrt((*eta_end) * (*eta_end) + xi_end[0] * xi_end[0] + xi_end[1] * xi_end[1] + xi_end[2] * xi_end[2]);
        for(int i = 0; i < 3; i++){
            xi_end[i] = xi_end[i] / sqrt((*eta_end) * (*eta_end) + xi_end[0] * xi_end[0] + xi_end[1] * xi_end[1] + xi_end[2] * xi_end[2]);
        }
        
        MatrixMulti_( 3, N, 1, J_g_v, q_dot, v_e);

		// ******************************************  计算末端走过的路程  ***************************************************************************
		(*locus) += sqrt(v_e[0] * v_e[0] + v_e[1] * v_e[1] + v_e[2] * v_e[2]) * delta_tau;
		//**************************************************************************************************************************************

        for(int i=0;i<3;i++){
            Pe[i] = Pe[i] + v_e[i] * delta_tau;
        }
        
        MatrixMulti_( 3, N, 1, J_bm_v, q_dot, v_b);
        for(int i=0;i<3;i++){
            r_b[i] = r_b[i] + v_b[i] * delta_tau;
        }

        quaternion2dc((*eta_b), xi_b[0], xi_b[1], xi_b[2], A_b);      
        links_transform(A_b, q, A_links_transform);       
        calc_r(r_b, A_b, A_links_transform, r);
        
        for(int i=0;i<3;i++){
		    b_tempt_N_minus_1[i] = b[(N-1)*3+i];		
		    for(int j=0;j<3;j++){
			    A_links_transform_N_minus_1[i*3+j] = A_links_transform[(N-1)*9+i*3+j];
            }
	    }
	    MatrixMulti_(3, 3, 1, A_links_transform_N_minus_1, b_tempt_N_minus_1, r_e_tempt);

        //  计算末端位置
	    for(int i=0;i<3;i++){
		    r_e[i] = r[(N-1)*3+i] + r_e_tempt[i]; 
        }
        
        // 计算 关节位置
	    calc_p(r, A_links_transform, p);
        
        //  Jm, Jbm, Jg
        for(int i=0;i<N;i++){
            //**************************************************************************************************************************************
		    MatrixExtract_( 1, 3*N, 1, 1, i*3+1, i*3+3, p, p_tempt);		
		    ScaleMatrix_( 1, 3, (-1), p_tempt, p_tempt);		
		    MatrixAdd_( 1, 3, r_e, p_tempt, r_e_minus_p_i);		
		    MatrixExtract_( 1, 3*3*N, 1, 1, i*9+1, i*9+9, A_links_transform, A_links_transform_tempt_i);
            //***************************************************************************************************************************************
		    
		    MatrixMulti_(3, 3, 1, A_links_transform_tempt_i, Ez, A_links_transform_tempt_i_multi_Ez);		    
		    cross(A_links_transform_tempt_i_multi_Ez, A_links_transform_tempt_i_multi_Ez_cross);		    
		    MatrixMulti_(3, 3, 1, A_links_transform_tempt_i_multi_Ez_cross, r_e_minus_p_i, Jm_v_tempt_i);
		    
		    for(int k = 0; k < 3; k++){
			    Jm_v_tempt[i * 3 + k] = Jm_v_tempt_i[k];
			    Jm_w_tempt[i * 3 + k] = A_links_transform_tempt_i_multi_Ez[k];
		    }
		    for(int jj = 0; jj < 6; jj++){
			    if(jj<3){
				    Jm_tempt[i * 6 + jj] = Jm_v_tempt[i * 3 + jj];
                }
			    else{
				    Jm_tempt[i * 6 + jj] = Jm_w_tempt[i * 3 + (jj - 3)];
                }
		    }
	    }
	
	    MatrixTranspose_(N, 3, Jm_v_tempt, Jm_v);
	    MatrixTranspose_(N, 3, Jm_w_tempt, Jm_w);
	    MatrixTranspose_(N, 6, Jm_tempt, Jm);
	    calc_J_bm(r, r_b, A_b, A_links_transform, p, J_bm_w, J_bm_v);

	    for(int j = 0; j < 6; j++){
		    for(int k = 0; k < N; k++){
				if(j < 3){
                    J_bm[j * N + k] = J_bm_v[j * N + k];
                }
			    else{
                    J_bm[j * N + k] = J_bm_w[(j-3) * N + k];
                }
		    }
        }
	    J_Base2EE(r_e, r_b, J_bE);
	    MatrixMulti_(6, 6, N, J_bE, J_bm, J_bE_multi_J_bm);
	    for(int i = 0; i < 6; i++){
            for(int j = 0; j < N; j++){
                J_g[i * N + j] = J_bE_multi_J_bm[i * N + j] + Jm[i * N + j];
            }
        }
		    
	    for(int i=0;i<3;i++){
		    for(int j=0;j<N;j++){
			    J_g_v[i*N+j] = J_g[i*N+j];
			    J_g_w[i*N+j] = J_g[(i+3)*N+j];
		    }
        } 

        MatrixMulti_(3, 3, nvrtx, A_b, base_center2vertex_temp, base_object_tempt);
        MatrixMulti_(3, 3, nvrtx, A_b, base_center2left_solar_vertex_temp, left_solar_object_tempt);
        MatrixMulti_(3, 3, nvrtx, A_b, base_center2right_solar_vertex_temp, right_solar_object_tempt);

        MatrixTranspose_(3, 8, base_object_tempt, base_object_tempt);
        MatrixTranspose_(3, 8, left_solar_object_tempt, left_solar_object_tempt);
        MatrixTranspose_(3, 8, right_solar_object_tempt, right_solar_object_tempt);

        MatrixExtract_( 1, 3*3*N, 1, 1, 0*9+1, 0*9+9, A_links_transform, A_links_transform_0);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_0, link_1_center2link_1_part_1_vertex_temp, link_1_part_1_object_tempt);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_0, link_1_center2link_1_part_2_vertex_temp, link_1_part_2_object_tempt);

        MatrixTranspose_(3, 8, link_1_part_1_object_tempt, link_1_part_1_object_tempt);
        MatrixTranspose_(3, 8, link_1_part_2_object_tempt, link_1_part_2_object_tempt);

        MatrixExtract_( 1, 3*3*N, 1, 1, 1*9+1, 1*9+9, A_links_transform, A_links_transform_1);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_1, link_2_center2link_2_part_1_vertex_temp, link_2_part_1_object_tempt);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_1, link_2_center2link_2_part_2_vertex_temp, link_2_part_2_object_tempt);

        MatrixTranspose_(3, 8, link_2_part_1_object_tempt, link_2_part_1_object_tempt);
        MatrixTranspose_(3, 8, link_2_part_2_object_tempt, link_2_part_2_object_tempt);

        MatrixExtract_( 1, 3*3*N, 1, 1, 4*9+1, 4*9+9, A_links_transform, A_links_transform_4);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_4, link_5_center2link_5_part_1_vertex_temp, link_5_part_1_object_tempt);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_4, link_5_center2link_5_part_2_vertex_temp, link_5_part_2_object_tempt);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_4, link_5_center2link_5_part_3_vertex_temp, link_5_part_3_object_tempt);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_4, link_5_center2link_5_part_4_vertex_temp, link_5_part_4_object_tempt);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_4, link_5_center2link_5_part_5_vertex_temp, link_5_part_5_object_tempt);

        MatrixTranspose_(3, 8, link_5_part_1_object_tempt, link_5_part_1_object_tempt);
        MatrixTranspose_(3, 8, link_5_part_2_object_tempt, link_5_part_2_object_tempt);
        MatrixTranspose_(3, 8, link_5_part_3_object_tempt, link_5_part_3_object_tempt);
        MatrixTranspose_(3, 8, link_5_part_4_object_tempt, link_5_part_4_object_tempt);
        MatrixTranspose_(3, 8, link_5_part_5_object_tempt, link_5_part_5_object_tempt);

        MatrixExtract_( 1, 3*3*N, 1, 1, 5*9+1, 5*9+9, A_links_transform, A_links_transform_5);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_5, link_6_center2link_6_part_1_vertex_temp, link_6_part_1_object_tempt);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_5, link_6_center2link_6_part_2_vertex_temp, link_6_part_2_object_tempt);

        MatrixTranspose_(3, 8, link_6_part_1_object_tempt, link_6_part_1_object_tempt);
        MatrixTranspose_(3, 8, link_6_part_2_object_tempt, link_6_part_2_object_tempt);

        MatrixExtract_( 1, 3*3*N, 1, 1, 6*9+1, 6*9+9, A_links_transform, A_links_transform_6);
        MatrixMulti_(3, 3, nvrtx, A_links_transform_6, link_7_center2link_7_vertex_temp, link_7_object_tempt);

        MatrixTranspose_(3, 8, link_7_object_tempt, link_7_object_tempt);

        for(int i = 0; i < nvrtx; i++){
            for(int j = 0; j < 3; j++){
                base_object[i][j]          = base_object_tempt[i*3 + j];
                left_solar_object[i][j]    = left_solar_object_tempt[i*3 + j];
                right_solar_object[i][j]   = right_solar_object_tempt[i*3 + j];
                link_1_part_1_object[i][j] = link_1_part_1_object_tempt[i*3 + j];
                link_1_part_2_object[i][j] = link_1_part_2_object_tempt[i*3 + j];
                link_2_part_1_object[i][j] = link_2_part_1_object_tempt[i*3 + j];
                link_2_part_2_object[i][j] = link_2_part_2_object_tempt[i*3 + j];
                link_5_part_1_object[i][j] = link_5_part_1_object_tempt[i*3 + j];
                link_5_part_2_object[i][j] = link_5_part_2_object_tempt[i*3 + j];
                link_5_part_3_object[i][j] = link_5_part_3_object_tempt[i*3 + j];
                link_5_part_4_object[i][j] = link_5_part_4_object_tempt[i*3 + j];
                link_5_part_5_object[i][j] = link_5_part_5_object_tempt[i*3 + j];
                link_6_part_1_object[i][j] = link_6_part_1_object_tempt[i*3 + j];
                link_6_part_2_object[i][j] = link_6_part_2_object_tempt[i*3 + j];
                link_7_object[i][j]        = link_7_object_tempt[i*3 + j];
            }
        }

        // base_object = (r_b + base_object_tempt).T
        for(int j = 0; j < 3; j++){
            for(int i = 0; i < nvrtx; i++){
                base_object[i][j] = base_object[i][j] + r_b[j];
                left_solar_object[i][j] = left_solar_object[i][j] + r_b[j];
                right_solar_object[i][j] = right_solar_object[i][j] + r_b[j];
                link_1_part_1_object[i][j] = link_1_part_1_object[i][j] + r[0 * 3 + j];
                link_1_part_2_object[i][j] = link_1_part_2_object[i][j] + r[0 * 3 + j];
                link_2_part_1_object[i][j] = link_2_part_1_object[i][j] + r[1 * 3 + j];
                link_2_part_2_object[i][j] = link_2_part_2_object[i][j] + r[1 * 3 + j];
                link_5_part_1_object[i][j] = link_5_part_1_object[i][j] + r[4 * 3 + j];
                link_5_part_2_object[i][j] = link_5_part_2_object[i][j] + r[4 * 3 + j];
                link_5_part_3_object[i][j] = link_5_part_3_object[i][j] + r[4 * 3 + j];
                link_5_part_4_object[i][j] = link_5_part_4_object[i][j] + r[4 * 3 + j];
                link_5_part_5_object[i][j] = link_5_part_5_object[i][j] + r[4 * 3 + j];
                link_6_part_1_object[i][j] = link_6_part_1_object[i][j] + r[5 * 3 + j];
                link_6_part_2_object[i][j] = link_6_part_2_object[i][j] + r[5 * 3 + j];
                link_7_object[i][j] = link_7_object[i][j] + r[6 * 3 + j];
            }
        }

        bd_base_object.coord = base_object;
        bd_base_object.numpoints = nvrtx;
        bd_left_solar_object.coord = left_solar_object;
        bd_left_solar_object.numpoints = nvrtx;
        bd_right_solar_object.coord = right_solar_object;
        bd_right_solar_object.numpoints = nvrtx;

        bd_link_1_part_1_object.coord = link_1_part_1_object;
        bd_link_1_part_1_object.numpoints = nvrtx;
        bd_link_1_part_2_object.coord = link_1_part_2_object;
        bd_link_1_part_2_object.numpoints = nvrtx;

        bd_link_2_part_1_object.coord = link_2_part_1_object;
        bd_link_2_part_1_object.numpoints = nvrtx;
        bd_link_2_part_2_object.coord = link_2_part_2_object;
        bd_link_2_part_2_object.numpoints = nvrtx;

        bd_link_5_part_1_object.coord = link_5_part_1_object;
        bd_link_5_part_1_object.numpoints = nvrtx;
        bd_link_5_part_2_object.coord = link_5_part_2_object;
        bd_link_5_part_2_object.numpoints = nvrtx;
        bd_link_5_part_3_object.coord = link_5_part_3_object;
        bd_link_5_part_3_object.numpoints = nvrtx;
        bd_link_5_part_4_object.coord = link_5_part_4_object;
        bd_link_5_part_4_object.numpoints = nvrtx;
        bd_link_5_part_5_object.coord = link_5_part_5_object;
        bd_link_5_part_5_object.numpoints = nvrtx;

        bd_link_6_part_1_object.coord = link_6_part_1_object;
        bd_link_6_part_1_object.numpoints = nvrtx;
        bd_link_6_part_2_object.coord = link_6_part_2_object;
        bd_link_6_part_2_object.numpoints = nvrtx;

        bd_link_7_object.coord = link_7_object;
        bd_link_7_object.numpoints = nvrtx;

        distance_base2link_5_part_1 = compute_minimum_distance(bd_base_object, bd_link_5_part_1_object, &s);
        distance_base2link_5_part_2 = compute_minimum_distance(bd_base_object, bd_link_5_part_2_object, &s);
        distance_base2link_5_part_3 = compute_minimum_distance(bd_base_object, bd_link_5_part_3_object, &s);
        distance_base2link_5_part_4 = compute_minimum_distance(bd_base_object, bd_link_5_part_4_object, &s);
        distance_base2link_5_part_5 = compute_minimum_distance(bd_base_object, bd_link_5_part_5_object, &s);
        distance_base2link_6_part_1 = compute_minimum_distance(bd_base_object, bd_link_6_part_1_object, &s);
        distance_base2link_6_part_2 = compute_minimum_distance(bd_base_object, bd_link_6_part_2_object, &s);
        distance_base2link_7        = compute_minimum_distance(bd_base_object, bd_link_7_object, &s);

        distance_left2link_5_part_2 = compute_minimum_distance(bd_left_solar_object, bd_link_5_part_2_object, &s);
        distance_left2link_5_part_3 = compute_minimum_distance(bd_left_solar_object, bd_link_5_part_3_object, &s);
        distance_left2link_5_part_4 = compute_minimum_distance(bd_left_solar_object, bd_link_5_part_4_object, &s);
        distance_left2link_5_part_5 = compute_minimum_distance(bd_left_solar_object, bd_link_5_part_5_object, &s);
        distance_left2link_6_part_1 = compute_minimum_distance(bd_left_solar_object, bd_link_6_part_1_object, &s);
        distance_left2link_6_part_2 = compute_minimum_distance(bd_left_solar_object, bd_link_6_part_2_object, &s);
        distance_left2link_7        = compute_minimum_distance(bd_left_solar_object, bd_link_7_object, &s);

        distance_right2link_5_part_2 = compute_minimum_distance(bd_right_solar_object, bd_link_5_part_2_object, &s);
        distance_right2link_5_part_3 = compute_minimum_distance(bd_right_solar_object, bd_link_5_part_3_object, &s);
        distance_right2link_5_part_4 = compute_minimum_distance(bd_right_solar_object, bd_link_5_part_4_object, &s);
        distance_right2link_5_part_5 = compute_minimum_distance(bd_right_solar_object, bd_link_5_part_5_object, &s);
        distance_right2link_6_part_1 = compute_minimum_distance(bd_right_solar_object, bd_link_6_part_1_object, &s);
        distance_right2link_6_part_2 = compute_minimum_distance(bd_right_solar_object, bd_link_6_part_2_object, &s);
        distance_right2link_7        = compute_minimum_distance(bd_right_solar_object, bd_link_7_object, &s);

        distance_link_1_part_1_link_5_part_2 = compute_minimum_distance(bd_link_1_part_1_object, bd_link_5_part_2_object, &s);
        distance_link_1_part_1_link_5_part_3 = compute_minimum_distance(bd_link_1_part_1_object, bd_link_5_part_3_object, &s);
        distance_link_1_part_1_link_5_part_4 = compute_minimum_distance(bd_link_1_part_1_object, bd_link_5_part_4_object, &s);
        distance_link_1_part_1_link_5_part_5 = compute_minimum_distance(bd_link_1_part_1_object, bd_link_5_part_5_object, &s);
        distance_link_1_part_1_link_6_part_1 = compute_minimum_distance(bd_link_1_part_1_object, bd_link_6_part_1_object, &s);
        distance_link_1_part_1_link_6_part_2 = compute_minimum_distance(bd_link_1_part_1_object, bd_link_6_part_2_object, &s);
        distance_link_1_part_1_link_7        = compute_minimum_distance(bd_link_1_part_1_object, bd_link_7_object, &s);

        distance_link_1_part_2_link_5_part_2 = compute_minimum_distance(bd_link_1_part_2_object, bd_link_5_part_2_object, &s);
        distance_link_1_part_2_link_5_part_3 = compute_minimum_distance(bd_link_1_part_2_object, bd_link_5_part_3_object, &s);
        distance_link_1_part_2_link_5_part_4 = compute_minimum_distance(bd_link_1_part_2_object, bd_link_5_part_4_object, &s);
        distance_link_1_part_2_link_5_part_5 = compute_minimum_distance(bd_link_1_part_2_object, bd_link_5_part_5_object, &s);
        distance_link_1_part_2_link_6_part_1 = compute_minimum_distance(bd_link_1_part_2_object, bd_link_6_part_1_object, &s);
        distance_link_1_part_2_link_6_part_2 = compute_minimum_distance(bd_link_1_part_2_object, bd_link_6_part_2_object, &s);
        distance_link_1_part_2_link_7        = compute_minimum_distance(bd_link_1_part_2_object, bd_link_7_object, &s);

        distance_link_2_part_1_link_5_part_2 = compute_minimum_distance(bd_link_2_part_1_object, bd_link_5_part_2_object, &s);
        distance_link_2_part_1_link_5_part_3 = compute_minimum_distance(bd_link_2_part_1_object, bd_link_5_part_3_object, &s);
        distance_link_2_part_1_link_5_part_4 = compute_minimum_distance(bd_link_2_part_1_object, bd_link_5_part_4_object, &s);
        distance_link_2_part_1_link_5_part_5 = compute_minimum_distance(bd_link_2_part_1_object, bd_link_5_part_5_object, &s);
        distance_link_2_part_1_link_6_part_1 = compute_minimum_distance(bd_link_2_part_1_object, bd_link_6_part_1_object, &s);
        distance_link_2_part_1_link_6_part_2 = compute_minimum_distance(bd_link_2_part_1_object, bd_link_6_part_2_object, &s);
        distance_link_2_part_1_link_7        = compute_minimum_distance(bd_link_2_part_1_object, bd_link_7_object, &s);

        distance_link_2_part_2_link_5_part_2 = compute_minimum_distance(bd_link_2_part_2_object, bd_link_5_part_2_object, &s);
        distance_link_2_part_2_link_5_part_3 = compute_minimum_distance(bd_link_2_part_2_object, bd_link_5_part_3_object, &s);
        distance_link_2_part_2_link_5_part_4 = compute_minimum_distance(bd_link_2_part_2_object, bd_link_5_part_4_object, &s);
        distance_link_2_part_2_link_5_part_5 = compute_minimum_distance(bd_link_2_part_2_object, bd_link_5_part_5_object, &s);
        distance_link_2_part_2_link_6_part_1 = compute_minimum_distance(bd_link_2_part_2_object, bd_link_6_part_1_object, &s);
        distance_link_2_part_2_link_6_part_2 = compute_minimum_distance(bd_link_2_part_2_object, bd_link_6_part_2_object, &s);
        distance_link_2_part_2_link_7        = compute_minimum_distance(bd_link_2_part_2_object, bd_link_7_object, &s);

        distance_base2link_567[0] = distance_base2link_5_part_1; distance_base2link_567[1] = distance_base2link_5_part_2; 
        distance_base2link_567[2] = distance_base2link_5_part_3; distance_base2link_567[3] = distance_base2link_5_part_4; 
        distance_base2link_567[4] = distance_base2link_5_part_5; distance_base2link_567[5] = distance_base2link_6_part_1; 
        distance_base2link_567[6] = distance_base2link_6_part_2; distance_base2link_567[7] = distance_base2link_7; 
        for(int ii = 0; ii < 8; ii++){
            if(distance_base2link_567[ii] > distance_limit){
                collision_test_base_link567 += 1; 
            }
        }
        distance_left2link_567[0] = distance_left2link_5_part_2; distance_left2link_567[1] = distance_left2link_5_part_3; 
        distance_left2link_567[2] = distance_left2link_5_part_4; distance_left2link_567[3] = distance_left2link_5_part_5; 
        distance_left2link_567[4] = distance_left2link_6_part_1; distance_left2link_567[5] = distance_left2link_6_part_2;
        distance_left2link_567[6] = distance_left2link_7;  
        for(int ii = 0; ii < 7; ii++){
            if(distance_left2link_567[ii] > distance_limit){
                collision_test_left_link567 += 1; 
            }
        }
        distance_right2link_567[0] = distance_right2link_5_part_2; distance_right2link_567[1] = distance_right2link_5_part_3;
        distance_right2link_567[2] = distance_right2link_5_part_4; distance_right2link_567[3] = distance_right2link_5_part_5;
        distance_right2link_567[4] = distance_right2link_6_part_1; distance_right2link_567[5] = distance_right2link_6_part_2;
        distance_right2link_567[6] = distance_right2link_7;
        for(int ii = 0; ii < 7; ii++){
            if(distance_right2link_567[ii] > distance_limit){
                collision_test_right_link567 += 1; 
            }
        }
        distance_link_12link_567[0]  = distance_link_1_part_1_link_5_part_2; distance_link_12link_567[1]  = distance_link_1_part_1_link_5_part_3; 
        distance_link_12link_567[2]  = distance_link_1_part_1_link_5_part_4; distance_link_12link_567[3]  = distance_link_1_part_1_link_5_part_5; 
        distance_link_12link_567[4]  = distance_link_1_part_1_link_6_part_1; distance_link_12link_567[5]  = distance_link_1_part_1_link_6_part_2; 
        distance_link_12link_567[6]  = distance_link_1_part_1_link_7;        distance_link_12link_567[7]  = distance_link_1_part_2_link_5_part_2; 
        distance_link_12link_567[8]  = distance_link_1_part_2_link_5_part_3; distance_link_12link_567[9]  = distance_link_1_part_2_link_5_part_4; 
        distance_link_12link_567[10] = distance_link_1_part_2_link_5_part_5; distance_link_12link_567[11] = distance_link_1_part_2_link_6_part_1; 
        distance_link_12link_567[12] = distance_link_1_part_2_link_6_part_2; distance_link_12link_567[13] = distance_link_1_part_2_link_7; 
        for(int ii = 0; ii < 14; ii++){
            if(distance_link_12link_567[ii] > distance_limit){
                collision_test_link_1_link567 += 1; 
            }
        }
        distance_link_22link_567[0]  = distance_link_2_part_1_link_5_part_2; distance_link_22link_567[1]  = distance_link_2_part_1_link_5_part_3; 
        distance_link_22link_567[2]  = distance_link_2_part_1_link_5_part_4; distance_link_22link_567[3]  = distance_link_2_part_1_link_5_part_5; 
        distance_link_22link_567[4]  = distance_link_2_part_1_link_6_part_1; distance_link_22link_567[5]  = distance_link_2_part_1_link_6_part_2; 
        distance_link_22link_567[6]  = distance_link_2_part_1_link_7;        distance_link_22link_567[7]  = distance_link_2_part_2_link_5_part_2; 
        distance_link_22link_567[8]  = distance_link_2_part_2_link_5_part_3; distance_link_22link_567[9]  = distance_link_2_part_2_link_5_part_4; 
        distance_link_22link_567[10] = distance_link_2_part_2_link_5_part_5; distance_link_22link_567[11] = distance_link_2_part_2_link_6_part_1; 
        distance_link_22link_567[12] = distance_link_2_part_2_link_6_part_2; distance_link_22link_567[13] = distance_link_2_part_2_link_7; 
        for(int ii = 0; ii < 14; ii++){
            if(distance_link_22link_567[ii] > distance_limit){
                collision_test_link_2_link567 += 1; 
            }
        }

        collision_test = collision_test + collision_test_base_link567 + collision_test_left_link567 + collision_test_right_link567 + collision_test_link_1_link567 + collision_test_link_2_link567;

    }

    //*************************************************************************************************************************************
    (*delta_xi_b_distrb_max) = delta_xi_base_norm_max;

    (*T_min) = T_min_temp;

    (*collision) = collision_test;    //  total collision times during whole motion process

	// ****************************************************    计算最终时刻的可操作度   ****************************************************
	
    double manipl_temp_2 = 0;
	MatrixDiagExpand( A_b, 3, 3, A_b_expand);
	MatrixTranspose_(6, 6, A_b_expand, A_b_expand_transpose);
	MatrixMulti_(6, 6, N, A_b_expand_transpose, J_g, A_b_expand_transpose_multi_J_g);
	MatrixTranspose_(6, N, A_b_expand_transpose_multi_J_g, A_b_expand_transpose_multi_J_g_transpose);
	MatrixMulti_(6, N, 6, A_b_expand_transpose_multi_J_g, A_b_expand_transpose_multi_J_g_transpose, manipl_temp);
	manipl_temp_2 = calc_determinantOfMatrix(manipl_temp, 6*6, 6, 6, 6);

	*manipl = sqrt(manipl_temp_2);




    delete[] base_object_tempt;
    delete[] left_solar_object_tempt;
    delete[] right_solar_object_tempt;
    delete[] link_1_part_1_object_tempt;
    delete[] link_1_part_2_object_tempt;
    delete[] link_2_part_1_object_tempt;
    delete[] link_2_part_2_object_tempt;
    delete[] link_5_part_1_object_tempt;
    delete[] link_5_part_2_object_tempt;
    delete[] link_5_part_3_object_tempt;
    delete[] link_5_part_4_object_tempt;
    delete[] link_5_part_5_object_tempt;
    delete[] link_6_part_1_object_tempt;
    delete[] link_6_part_2_object_tempt;
    delete[] link_7_object_tempt;
    delete[] A_links_transform_0;
    delete[] A_links_transform_1;
    delete[] A_links_transform_4;
    delete[] A_links_transform_5;
    delete[] A_links_transform_6;

    delete[] base_center2vertex_temp;
    delete[] base_center2left_solar_vertex_temp;
    delete[] base_center2right_solar_vertex_temp;
    delete[] link_1_center2link_1_part_1_vertex_temp;
    delete[] link_1_center2link_1_part_2_vertex_temp;
    delete[] link_2_center2link_2_part_1_vertex_temp;
    delete[] link_2_center2link_2_part_2_vertex_temp;
    delete[] link_5_center2link_5_part_1_vertex_temp;
    delete[] link_5_center2link_5_part_2_vertex_temp;
    delete[] link_5_center2link_5_part_3_vertex_temp;
    delete[] link_5_center2link_5_part_4_vertex_temp;
    delete[] link_5_center2link_5_part_5_vertex_temp;
    delete[] link_6_center2link_6_part_1_vertex_temp;
    delete[] link_6_center2link_6_part_2_vertex_temp;
    delete[] link_7_center2link_7_vertex_temp;

    delete[] distance_base2link_567;
    delete[] distance_left2link_567;
    delete[] distance_right2link_567;
    delete[] distance_link_12link_567;
    delete[] distance_link_22link_567;



    for (int i = 0; i < nvrtx; i++){
        delete[] (bd_base_object.coord[i]);
        delete[] (bd_left_solar_object.coord[i]);
        delete[] (bd_right_solar_object.coord[i]);
        delete[] (bd_link_1_part_1_object.coord[i]);
        delete[] (bd_link_1_part_2_object.coord[i]);
        delete[] (bd_link_2_part_1_object.coord[i]);
        delete[] (bd_link_2_part_2_object.coord[i]);
        delete[] (bd_link_5_part_1_object.coord[i]);
        delete[] (bd_link_5_part_2_object.coord[i]);
        delete[] (bd_link_5_part_3_object.coord[i]);
        delete[] (bd_link_5_part_4_object.coord[i]);
        delete[] (bd_link_5_part_5_object.coord[i]);
        delete[] (bd_link_6_part_1_object.coord[i]);
        delete[] (bd_link_6_part_2_object.coord[i]);
        delete[] (bd_link_7_object.coord[i]);
    }
    delete[] (bd_base_object.coord);    
    delete[] (bd_left_solar_object.coord);
    delete[] (bd_right_solar_object.coord);    
    delete[] (bd_link_1_part_1_object.coord);
    delete[] (bd_link_1_part_2_object.coord);
    delete[] (bd_link_2_part_1_object.coord);    
    delete[] (bd_link_2_part_2_object.coord);
    delete[] (bd_link_5_part_1_object.coord);    
    delete[] (bd_link_5_part_2_object.coord);
    delete[] (bd_link_5_part_3_object.coord);
    delete[] (bd_link_5_part_4_object.coord);
    delete[] (bd_link_5_part_5_object.coord);
    delete[] (bd_link_6_part_1_object.coord);
    delete[] (bd_link_6_part_2_object.coord);
    delete[] (bd_link_7_object.coord);


	delete[] q ;
    delete[] q_dot ;
    delete[] q_ddot ;
    delete[] quaternion_base_initial ;
    delete[] A_b ;
    delete[] xi_b_tempt ;
    delete[] xi_end_tempt ;
    delete[] A_links_transform ;
    delete[] A_b_multi_b_b ;
    delete[] r_b_tempt_1 ;
    delete[] r_b_tempt_2 ;
    delete[] r_b_tempt_2_temp ;
    delete[] a_tempt_i ;
    delete[] r_b_tempt_3 ;
    delete[] r_b ;
    delete[] a_tempt_k ;
    delete[] b_tempt_k ;
    delete[] l_tempt_k ;    
    delete[] A_links_transform_tempt_i ;
    delete[] A_links_transform_tempt_i_multi_a_tempt_i;
    delete[] A_links_transform_tempt_k ;
    delete[] A_links_transform_tempt_k_multi_l_tempt_k ;
    delete[] r_b_tempt_3_tempt ;
    delete[] r ;
    delete[] r_e ;
    delete[] Pe_initial ;
    delete[] A_links_transform_N_minus_1;
    delete[] r_e_tempt;
    delete[] b_tempt_N_minus_1 ;
    delete[] p ;
    delete[] A_links_transform_tempt_i_multi_Ez ;
    delete[] A_links_transform_tempt_i_multi_Ez_cross;
    delete[] Jm_v ;
    delete[] Jm_v_tempt ;
    delete[] Jm_v_tempt_i ;
    delete[] Jm_w ;
    delete[] Jm_w_tempt ;
    delete[] p_tempt ;
    delete[] r_e_minus_p_i ;
    delete[] Jm ;
    delete[] Jm_tempt;
    delete[] J_bm_w ;
    delete[] J_bm_v ;
    delete[] J_bm ;	
    delete[] J_bE ;
    delete[] J_g ;
    delete[] J_bE_multi_J_bm ;
    delete[] J_g_v ;
    delete[] J_g_w ;    
    delete[] quaternion_end ;
    delete[] eta_b_dot_tempt_1 ;
    delete eta_b_dot_tempt_2 ;  
    delete[] xi_b_dot ;
    delete[] cross_xi_b ;  
    delete[] xi_b_dot_tempt_1 ;
    delete[] xi_b_dot_tempt_2 ;
    delete[] xi_b_dot_tempt_3 ; 
    delete[] eta_end_dot_tempt_1 ;
    delete eta_end_dot_tempt_2 ;
    delete[] xi_end_dot ;
    delete[] cross_xi_end; 
    delete[] xi_end_dot_tempt_1 ;
    delete[] xi_end_dot_tempt_2 ;
    delete[] xi_end_dot_tempt_3 ;
    delete[] v_e ;
    delete[] v_b ;
    delete[] A_b_expand ;
	delete[] A_b_expand_transpose ;
	delete[] A_b_expand_transpose_multi_J_g ;
	delete[] A_b_expand_transpose_multi_J_g_transpose ;
	delete[] manipl_temp ;
    delete[] xi_b_initial ;
    delete[] delta_xi_base ;
    delete delta_eta_base_tempt;
    delete delta_eta_base ;
    delete[] delta_xi_base_tempt;


}
























