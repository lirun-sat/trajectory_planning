#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

int main()
{
	// double para[] = {1.33997, 0.98804, -3.811, -1.77003, -2.93483, -0.321918, -1.28835};
	// double para[] = {1.82072 , 0.587827 , 1.0189 , 4.33098 , -3.20266 , 1.06863 , 5.27103};
	// double para[] = {1.23848 , 0.400183 , -0.838996 , 2.94077 , 1.44408 , -2.0944 , -1.96748};
    // double para[] = {1.143, 0.322779, 6.28319, 3.1469, 3.81277, 2.0944, 1.05599};
	// double para[] = {0.870048, 0.340917, 2.57864, 2.91813, -1.85878, -1.91102, -2.1036};
    double para[] = {0.826199, 0.384593, 2.0656, 2.86412, 1.74332, 1.9403, 0.964951};

    double* eta_end = new double;
    double* xi_end = new double[3];
    double* Pe = new double[3];
    double* eta_b = new double;
    double* xi_b = new double[3];
    double* quaternion_end_desired = new double[4];
    double* eta_base_initial = new double;
    double* xi_base_initial = new double[3];
    double* eta_end_desired = new double;
    double* xi_end_desired = new double[3];
    double* delta_eta_base;
    double* delta_eta_end;
    double* delta_xi_base;
    double* delta_xi_end;
    double* delta_Pe_end;
    double* rpy_end = new double[3];
    double* rpy_base = new double[3];
    double* delta_euler_end = new double[3];
    double* delta_euler_base = new double[3];

    double* p_e_initial = new double[3];
    double* locus = new double;
    double* delta_xi_b_distrb_max = new double;
    double* manipl = new double;
    double* T_min = new double;
    
    delta_eta_base = new double;
    delta_eta_end = new double;
    delta_xi_base = new double[3];
    delta_xi_end = new double[3];
    delta_Pe_end = new double[3];
    
    double K_a = 1 / 0.0002;  // 1 degree error tolerance
	double K_p = 1 / 0.002;  //  0.002 meter error tolerance for end-effector
    double K_b = 1 / 0.0008;  // attitude error tolerance for base is 5 degree
    double K_M = 1;
    double K_t = 1 / 100;  // 1 / 10;  // 1;  // 1 / 10;  //  max allowed motion time is set to 100 seconds, tolerance is 10 seconds.

	double cost_func = 0;
    double delta_xi_end_mod_temp = 0;   
    double delta_Pe_end_mod_temp = 0;
    double delta_xi_base_mod_temp = 0;
    
	forward_kin_2(para, eta_end,  xi_end, Pe, eta_b, xi_b, p_e_initial, locus, delta_xi_b_distrb_max, manipl, T_min);

    // ********************************************************  straight_line_locus  ************************************************************
    // double delta_p_e = 0;
    // double straight_line_locus = 0;
    // for(int i = 0; i < 3; i++)
    // {
    //     delta_p_e = p_e_initial[i] - Pe[i];
    //     straight_line_locus += delta_p_e * delta_p_e;
    // }
    // straight_line_locus = sqrt(straight_line_locus);
// **********************************************************************************************************************************    
    cout << "delta_tau: " << delta_tau << endl;
    
    cout << "p_e_initial: " << endl;
    cout << p_e_initial[0] << "  " << p_e_initial[1] << "  " << p_e_initial[2] << endl;

    cout << "Pe: " << endl;
    cout << Pe[0] << "  " << Pe[1] << "  " << Pe[2] << endl;

    cout << "quaternion_end:" << endl;
    cout << *eta_end << "  " << xi_end[0] << "  " << xi_end[1] << "  " << xi_end[2] << endl;

    zyx2quaternion( RPY_END_DESIRED[2],  RPY_END_DESIRED[1],  RPY_END_DESIRED[0],  quaternion_end_desired);
    
    *eta_end_desired = quaternion_end_desired[0];

    for(int i = 0; i < 3; i++)
    {
        xi_end_desired[i] = quaternion_end_desired[i+1];
    }
    
    cout << "quaternion_end_desired:" << endl;
    cout << *eta_end_desired << "  " << xi_end_desired[0] << "  " << xi_end_desired[1] << "  " << xi_end_desired[2] << endl;
    
    delta_var(Pe_DESIRED, eta_end_desired, xi_end_desired, Pe, *eta_end, xi_end, *eta_b, xi_b, delta_eta_end, delta_xi_end, delta_Pe_end, delta_eta_base, delta_xi_base);    
    
    for(int i = 0; i < 3; i++)
    {
        delta_xi_end_mod_temp += delta_xi_end[i] * delta_xi_end[i];
        delta_Pe_end_mod_temp += delta_Pe_end[i] * delta_Pe_end[i];
        delta_xi_base_mod_temp += delta_xi_base[i] * delta_xi_base[i];
    }

    delta_xi_end_mod_temp = sqrt(delta_xi_end_mod_temp);
    delta_Pe_end_mod_temp = sqrt(delta_Pe_end_mod_temp);
    delta_xi_base_mod_temp = sqrt(delta_xi_base_mod_temp);

    // ****************************************  RPY error of end-effector,  Pe error  ********************************************************************************
    cost_func = K_a * delta_xi_end_mod_temp + K_p * delta_Pe_end_mod_temp + K_b * (*delta_xi_b_distrb_max) + K_M * (1 / (*manipl)) + K_t * (*T_min);


    cout << "cost_func: " << cost_func << endl; 
    cout << "manipl: " << *manipl << endl; 
    cout << "T_min: " << *T_min << endl;
    cout << "eta_b: " << *eta_b << endl; 
    cout << "xi_b[0]: " << xi_b[0] << endl;
    cout << "xi_b[1]: " << xi_b[1] << endl; 
    cout << "xi_b[2]: " << xi_b[2] << endl;


    

    delete[] delta_eta_base;
    delete[] delta_eta_end;
    delete[] delta_xi_base;
    delete[] delta_xi_end;
    delete[] delta_Pe_end;
    delete[] xi_end_desired;
    delete[] eta_end_desired;
    delete[] quaternion_end_desired;
    delete[] delta_euler_base;
    delete[] eta_base_initial;
    delete[] xi_base_initial;
    delete[] delta_euler_end;
    delete[] xi_b;
    delete[] eta_b;
    delete[] Pe;
    delete[] xi_end;
    delete[] eta_end;
    delete[] rpy_end;
    delete[] rpy_base;
    delete[] p_e_initial;
    delete locus;
    delete delta_xi_b_distrb_max;
    delete manipl;
    delete T_min;
	
}
