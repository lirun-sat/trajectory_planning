#ifndef CALC_FORWARD_KINEMATICS
#define CALC_FORWARD_KINEMATICS


void A_b(double alpha_base, double beta_base, double gamma_base, double* A_base);

void links_transform(double* A_base, double* q, double* A_links_transform);

void calc_r(double* r_b, double* A_b, double* A_links_transform, double* r);

void calc_p(double* r, double* A_links_transform, double* p);

void calc_J_bm(double* r, double* r_b, double* A_b, double* A_links_transform, double* p, double* J_bm_w, double* J_bm_v);

void J_Base2EE(double* r_e, double* r_b, double* J_bE);

double calc_q_dot(double tau);
double calc_q_ddot(double tau);

int calc_binomial(int n, int k);

void calc_para_range(double* result_min, double* result_max);

void forward_kin(double* para, double* eta_end, double* xi_end, double* Pe, double* eta_b, double* xi_b);

void forward_kin_2(double* para, double* eta_end, double* xi_end, double* Pe, double* eta_b, double* xi_b, 
				   double* p_e_initial, double* locus, double* RPY_base_distrb_max, double* manipl, double* T_min);

void delta_var(double* Pe_desired, double* eta_end_desired, double* xi_end_desired, 
			   double* Pe, double eta_end, double* xi_end, double eta_b, double* xi_b,
			   double* delta_eta_end, double* delta_xi_end, double* delta_Pe_end, 
			   double* delta_eta_base, double* delta_xi_base);

int calc_binomial(int n, int k);



#endif //CALC_FORWARD_KINEMATICS









