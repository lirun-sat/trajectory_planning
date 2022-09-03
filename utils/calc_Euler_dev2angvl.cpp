#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

void calc_Euler_dev2angvl(double* Euler_angle, double* Euler_dev, double* angvl)
{
    /*
    N_Phi = np.array([
        [0, -sin(alpha_z), cos(alpha_z) * cos(beta_y)],
        [0,  cos(alpha_z), sin(alpha_z) * cos(beta_y)],
        [1,  0,            -sin(beta_y)]
    ])
    Euler_dev = np.array([
        [alpha_dot],
        [beta_dot],
        [gamma_dot]
    ])
    omega = N_Phi @ Euler_dev
    */
    
    double* N_phi;
    N_phi = new double[3*3];
    
    double alpha_z = Euler_angle[0];
    double beta_y  = Euler_angle[1];
    double gamma_x = Euler_angle[2];
    
    N_phi[0] = 0; N_phi[1] = -sin(alpha_z); N_phi[2] = cos(alpha_z) * cos(beta_y);
    N_phi[3] = 0; N_phi[4] = cos(alpha_z);  N_phi[5] = sin(alpha_z) * cos(beta_y);
    N_phi[6] = 1; N_phi[7] = 0;             N_phi[8] = -sin(beta_y);
    
    MatrixMulti_( 3,  3,  1,  N_phi,  Euler_dev,  angvl);
    
    delete[] N_phi;
    
}
