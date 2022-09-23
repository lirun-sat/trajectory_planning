#define _USE_MATH_DEFINES 
#include <cmath>


// double delta_tau = 0.01;
// double delta_tau = 0.005;
double delta_tau = 0.001;


double rpy_joints[] = {
      0,      0, M_PI/2,
-M_PI/2, M_PI/2,      0,
 M_PI/2,      0,   M_PI,
-M_PI/2,      0,      0,
 M_PI/2,      0,   M_PI,
 M_PI/2,      0,      0,
-M_PI/2,      0,      0};


/*
double q_INITIAL[] = {
0, 
0.78, 
1.57,
0.78, 
0, 
-1.57, 
0};
*/


double q_INITIAL[] = {
0, 
0, 
0,
0, 
0, 
0, 
0};

double RPY_END_INITIAL[] = {0, M_PI / 2, M_PI / 2};                                    
double RPY_BASE_INITIAL[] = {0, 0, 0};
double RPY_END_DESIRED[] = {-30 * M_PI / 180, -50 * M_PI /180, 45 * M_PI /180 };
double Pe_DESIRED[] = {1.5, -1.4, 2.2};

int N = 7;

double Ez[] = {0, 0, 1};
double eye[] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; 
double joint_angle_velocity_min_limit = -0.0873;
double joint_angle_velocity_max_limit = 0.0873;    // angle velocity 5 deg/s
double joint_angle_acceleration_min_limit = -0.00873;
double joint_angle_acceleration_max_limit = 0.00873;

double joint_angle_min_limit_degree[] = {-360, 0, -360, -360, -360, -120, -360};

double joint_angle_min_limit_rad[] = {joint_angle_min_limit_degree[0] * M_PI / 180,
                                      joint_angle_min_limit_degree[1] * M_PI / 180,
                                      joint_angle_min_limit_degree[2] * M_PI / 180,
                                      joint_angle_min_limit_degree[3] * M_PI / 180,
                                      joint_angle_min_limit_degree[4] * M_PI / 180,
                                      joint_angle_min_limit_degree[5] * M_PI / 180,
                                      joint_angle_min_limit_degree[6] * M_PI / 180};


double joint_angle_max_limit_degree[] = {360, 180, 360, 360, 360, 120, 360};

double joint_angle_max_limit_rad[] = {joint_angle_max_limit_degree[0] * M_PI / 180,
                                      joint_angle_max_limit_degree[1] * M_PI / 180,
                                      joint_angle_max_limit_degree[2] * M_PI / 180,
                                      joint_angle_max_limit_degree[3] * M_PI / 180,
                                      joint_angle_max_limit_degree[4] * M_PI / 180,
                                      joint_angle_max_limit_degree[5] * M_PI / 180,
                                      joint_angle_max_limit_degree[6] * M_PI / 180};


double m_b = 1474.1;
double b_b[] = {0, 0, 1.0864};
double m[] = {39.618, 13.427, 22.757, 13.427, 42.614, 13.427, 8.269};

double a[] = {
0, 0.013631, 0.13304,
0, 0.018673, 0.08,
0, 0.011017, 0.8422,
0, 0.018673, 0.08,
0, 0.113430, 1.1109,
0, 0.018673, 0.08,
0, 0, 0.10565};

double b[] = {
0, 0.136369, 0.056960,
0, 0.1-0.018673, 0,
0, 0.1-0.011017, 1.08-0.8422,
0, 0.1-0.018673, 0,
0, 0.08-0.11343, 1.26-1.1109,
0, 0.1-0.018673, 0,
0, 0, 0.16435};

double I_b_body[] = {
17388.34, 0, 0,
0, 1340.43, 0,
0, 0, 17981.26};

double I_links_body[] = {
0.30554, 0, 0,
0, 0.25403, 0.030757,
0, 0.030757, 0.15563,

0.043493, 0, 0,
0, 0.031896, 0,
0, 0, 0.029035,

2.69, 0, 0,
0, 2.68, 0,
0, 0, 0.06,

0.043493, 0, 0,
0, 0.031896, 0,
0, 0, 0.029035,

1.75, 0, 0,
0, 1.47, 0.29,
0, 0.29, 0.33,

0.043493, 0, 0,
0, 0.031896, 0,
0, 0, 0.029035,

0.04, 0, 0,
0, 0.04, 0,
0, 0, 0.01};


// double base_center[] = {0, 0, 1.0136};
// double link_1_center[] = {0, 0.013631, 0.13304};
// double link_2_center[] = {0, 0.018673, 0.08};
// double link_5_center[] = {0, 0.11343, 1.1109};
// double link_6_center[] = {0, 0.018673, 0.08};
// double link_7_center[] = {0, 0, 0.10565};


int nvrtx = 8;

double base_center2vertex[][3] = {
{-1.0000,  1.0000 ,  0.9864},                            
{1.0000 ,  1.0000 ,  0.9864},
{1.0000 ,  -1.0000,  0.9864},
{-1.0000,  -1.0000,  0.9864},
{-1.0000,  1.0000 , -1.0136},
{1.0000 ,  1.0000 , -1.0136},
{1.0000 ,  -1.0000, -1.0136},
{-1.0000,  -1.0000, -1.0136}};


double base_center2left_solar_vertex[][3] = {
{1.5000, 1.0000,   0.0364},
{1.5000 , 1.0000,  -0.0636},
{-1.5000, 1.0000,   0.0364},
{-1.5000, 1.0000,   -0.0636},
{1.5000 , 7.2000,   0.0364},
{1.5000 , 7.2000,   -0.0636},
{-1.5000, 7.2000,   0.0364},
{-1.5000, 7.2000,   -0.0636}};

double base_center2right_solar_vertex[][3] = {
{1.5000 , -7.2000 , 0.0364},
{1.5000 , -7.2000 , -0.0636},
{-1.5000,  -7.2000,  0.0364},
{-1.5000,  -7.2000,  -0.0636},
{1.5000 , -1.0000 , 0.0364},
{1.5000 , -1.0000 , -0.0636},
{-1.5000,  -1.0000,  0.0364},
{-1.5000,  -1.0000,  -0.0636}};

double link_1_center2link_1_part_1_vertex[][3] = {
{0.075000 , -0.088631 ,  -0.133040},
{0.075000 , 0.061369  ,  -0.133040},
{-0.075000,  -0.088631,  -0.133040},
{-0.075000,  0.061369 ,  -0.133040},
{0.075000 , -0.088631 ,   0.116960},
{0.075000 , 0.061369  ,   0.116960},
{-0.075000,  -0.088631,   0.116960},
{-0.075000,  0.061369 ,   0.116960}};

double link_1_center2link_1_part_2_vertex[][3] = {
{0.050000 ,  0.061369 , 0.006960},
{0.050000 ,  0.061369 , 0.106960},
{-0.050000,  0.061369 , 0.106960},
{-0.050000,  0.061369 , 0.006960},
{-0.050000,  0.136369 , 0.006960},
{-0.050000,  0.136369 , 0.106960},
{0.050000 ,  0.136369 , 0.106960},
{0.050000 ,  0.136369 , 0.006960}};

double link_2_center2link_2_part_1_vertex[][3] = {
{0.050000 , -0.068673 ,  -0.080000},
{0.050000 , 0.031327  ,  -0.080000},
{-0.050000,  0.031327 ,  -0.080000},
{-0.050000,  -0.068673,  -0.080000},
{0.050000 , -0.068673 ,   0.080000},
{0.050000 , 0.031327  ,   0.080000},
{-0.050000,  0.031327 ,   0.080000},
{-0.050000,  -0.068673,   0.080000}};

double link_2_center2link_2_part_2_vertex[][3] = {
{0.050000 , 0.031327 ,  -0.050000},
{0.050000 , 0.081327 , -0.050000},
{-0.050000,  0.081327,  -0.050000},
{-0.050000,  0.031327,  -0.050000},
{0.050000 , 0.031327 , 0.050000},
{0.050000 , 0.081327 , 0.050000},
{-0.050000,  0.081327,  0.050000},
{-0.050000,  0.031327,  0.050000}};

double link_5_center2link_5_part_1_vertex[][3] = {
{0.030000 , -0.143430 ,  -1.110900},
{0.030000 , -0.083430 , -1.110900},
{-0.030000,  -0.143430,  -1.110900},
{-0.030000,  -0.083430,  -1.110900},
{0.030000 , -0.143430 , -0.110900},
{0.030000 , -0.083430 , -0.110900},
{-0.030000,  -0.143430,  -0.110900},
{-0.030000,  -0.083430,  -0.110900}};

double link_5_center2link_5_part_2_vertex[][3] = {
{0.050000 , -0.163430 ,  -0.110900},
{0.050000 , -0.063430 , -0.110900},
{-0.050000,  -0.163430,  -0.110900},
{-0.050000,  -0.063430,  -0.110900},
{0.050000 , -0.163430 , 0.049100},
{0.050000 , -0.063430 , 0.049100},
{-0.050000,  -0.163430,  0.049100},
{-0.050000,  -0.063430,  0.049100}};

double link_5_center2link_5_part_3_vertex[][3] = {
{0.050000 , -0.063430 ,  -0.08090},
{-0.050000,  -0.063430,  -0.08090},
{0.050000 , -0.063430 ,  0.019100},
{-0.050000,  -0.063430,  0.019100},
{0.050000 , 0.146570  , -0.080900},
{-0.050000,  0.146570 , -0.080900},
{0.050000 , 0.146570  ,  0.019100},
{-0.050000,  0.146570 ,  0.019100}};

double link_5_center2link_5_part_4_vertex[][3] = {
{0.050000 , 0.016570 ,  0.019100},
{0.050000 , 0.116570 , 0.019100},
{-0.050000,  0.016570,  0.019100},
{-0.050000,  0.116570,  0.019100},
{0.050000 , 0.016570 , 0.229100},
{0.050000 , 0.116570 , 0.229100},
{-0.050000,  0.016570,  0.229100},
{-0.050000,  0.116570,  0.229100}};

double link_5_center2link_5_part_5_vertex[][3] = {
{0.050000 , 0.016570  ,  0.099100},
{-0.050000,  0.016570 ,  0.099100},
{-0.050000,  0.016570 ,  0.199100},
{0.050000 , 0.016570  ,  0.199100},
{0.050000 , -0.033430 ,  0.099100},
{-0.050000,  -0.033430,  0.099100},
{-0.050000,  -0.033430,  0.199100},
{0.050000 , -0.033430 , 0.199100}};

double link_6_center2link_6_part_1_vertex[][3] = {
{0.050000 , -0.068673 ,  -0.080000},
{-0.050000,  -0.068673,  -0.080000},
{-0.050000,  0.031327 ,  -0.080000},
{0.050000 , 0.031327  , -0.080000},
{0.050000 , -0.068673 ,  0.080000},
{-0.050000,  -0.068673,  0.080000},
{-0.050000,  0.031327 ,  0.080000},
{0.050000 , 0.031327  ,  0.080000}};

double link_6_center2link_6_part_2_vertex[][3] = {
{0.050000 , 0.031327 ,  -0.050000},
{-0.050000,  0.031327,  -0.050000},
{-0.050000,  0.031327,  0.050000},
{0.050000 , 0.031327 , 0.050000},
{0.050000 , 0.081327 , -0.050000},
{-0.050000,  0.081327,  -0.050000},
{-0.050000,  0.081327,  0.050000},
{0.050000 , 0.081327 , 0.050000}};

double link_7_center2link_7_vertex[][3] = {
{-0.040000,  -0.040000, -0.105650},
{-0.040000,  0.040000 , -0.105650},
{0.040000 , 0.040000  , -0.105650},
{0.040000 , -0.040000 , -0.105650},
{-0.040000,  -0.040000,  0.164350},
{-0.040000,  0.040000 ,  0.164350},
{0.040000 , 0.040000  ,  0.164350},
{0.040000 , -0.040000 ,  0.164350}};

