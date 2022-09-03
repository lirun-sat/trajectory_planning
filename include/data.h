extern int N;
extern double Ez[];
extern double eye[];

extern double joint_angle_min_limit_degree[];
extern double joint_angle_min_limit_rad[];
extern double joint_angle_max_limit_degree[];
extern double joint_angle_max_limit_rad[];

extern double joint_angle_velocity_min_limit;
extern double joint_angle_velocity_max_limit;
extern double joint_angle_acceleration_min_limit;
extern double joint_angle_acceleration_max_limit;

extern double rpy_joints[];
extern double m_b;
extern double b_b[];
extern double m[];
extern double a[];
extern double b[];
extern double I_b_body[];
extern double I_links_body[];
extern double q_INITIAL[];
extern double RPY_END_INITIAL[];
extern double delta_tau;
extern double Pe_DESIRED[];
extern double RPY_END_DESIRED[];
extern double RPY_BASE_INITIAL[];
//******************************************************************************************************************************
extern double base_vertex[][3];
extern double left_solar_vertex[][3];
extern double right_solar_vertex[][3];
extern double link_1_part_1_vertex[][3];
extern double link_1_part_2_vertex[][3];
extern double link_2_part_1_vertex[][3];
extern double link_2_part_2_vertex[][3];
extern double link_5_part_1_vertex[][3];
extern double link_5_part_2_vertex[][3];
extern double link_5_part_3_vertex[][3];
extern double link_5_part_4_vertex[][3];
extern double link_5_part_5_vertex[][3];
extern double link_6_part_1_vertex[][3];
extern double link_6_part_2_vertex[][3];
extern double link_7_vertex[][3];

extern double base_center[3];
extern double link_1_center[3];
extern double link_2_center[3];
extern double link_5_center[3];
extern double link_6_center[3];
extern double link_7_center[3];
//*****************************************************************************************************************************
extern double base_center2vertex[][3];
extern double base_center2left_solar_vertex[][3];
extern double base_center2right_solar_vertex[][3];

extern double link_1_center2link_1_part_1_vertex[][3];
extern double link_1_center2link_1_part_2_vertex[][3];

extern double link_2_center2link_2_part_1_vertex[][3];
extern double link_2_center2link_2_part_2_vertex[][3];

extern double link_5_center2link_5_part_1_vertex[][3];
extern double link_5_center2link_5_part_2_vertex[][3];
extern double link_5_center2link_5_part_3_vertex[][3];
extern double link_5_center2link_5_part_4_vertex[][3];
extern double link_5_center2link_5_part_5_vertex[][3];

extern double link_6_center2link_6_part_1_vertex[][3];
extern double link_6_center2link_6_part_2_vertex[][3];

extern double link_7_center2link_7_vertex[][3];

extern int nvrtx;


