#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void delta_var(double* Pe_desired, double* eta_end_desired, double* xi_end_desired, double* Pe, double eta_end, double* xi_end, double eta_b, double* xi_b,
                double* delta_eta_end, double* delta_xi_end, double* delta_Pe_end, double* delta_eta_base, double* delta_xi_base)
{
    //计算最终末端位姿与期望位姿差值,基座最终姿态与初始姿态之间的差值
    double* quaternion_base_initial = new double[4];
    double eta_b_initial;
    double* xi_b_initial = new double[3];
    double* delta_xi_base_tempt = new double[3];
    double* cross_xi_b = new double[3*3];
    double* delta_eta_base_tempt = new double;
    double* cross_xi_end = new double[3*3];
    double* delta_xi_end_tempt = new double[3];
    double* delta_eta_end_tempt = new double;

    zyx2quaternion(RPY_BASE_INITIAL[2], RPY_BASE_INITIAL[1], RPY_BASE_INITIAL[0], quaternion_base_initial);
    eta_b_initial = quaternion_base_initial[0];
    for (int i = 0; i < 3; i++)
    {
        xi_b_initial[i] = quaternion_base_initial[i + 1];
    }
        
    
    cross(xi_b, cross_xi_b);
    MatrixMulti_(3, 3, 1, cross_xi_b, xi_b_initial, delta_xi_base_tempt);
    for(int i = 0; i < 3; i++)
    {
        delta_xi_base_tempt[i] = (-1) * delta_xi_base_tempt[i] - eta_b_initial * xi_b[i];     
        delta_xi_base_tempt[i] = delta_xi_base_tempt[i] + eta_b * xi_b_initial[i];
    }
    
    for (int i = 0; i < 3; i++)
    {
        delta_xi_base[i] = delta_xi_base_tempt[i];
    }
        
    MatrixMulti_(1, 3, 1, xi_b, xi_b_initial, delta_eta_base_tempt);
    (*delta_eta_base_tempt) = (*delta_eta_base_tempt) + eta_b * eta_b_initial;
    *delta_eta_base = (*delta_eta_base_tempt);

    cross(xi_end, cross_xi_end);
    MatrixMulti_(3, 3, 1, cross_xi_end, xi_end_desired, delta_xi_end_tempt);
    for(int i = 0; i < 3; i++)
    {
        delta_xi_end_tempt[i] = (-1) * delta_xi_end_tempt[i] - (*eta_end_desired) * xi_end[i];        
        delta_xi_end_tempt[i] = delta_xi_end_tempt[i] + eta_end * xi_end_desired[i];
    }
    
    for (int i = 0; i < 3; i++)
    {
        delta_xi_end[i] = delta_xi_end_tempt[i];
    }
    MatrixMulti_(1, 3, 1, xi_end, xi_end_desired, delta_eta_end_tempt);
    (*delta_eta_end_tempt) = (*delta_eta_end_tempt) + eta_end * (*eta_end_desired);    
    (*delta_eta_end) = (*delta_eta_end_tempt);
    
    for (int i = 0; i < 3; i++)
    {
        delta_Pe_end[i] = Pe_desired[i] - Pe[i];
    }


    delete[] quaternion_base_initial;
    delete[] xi_b_initial;
    delete[] delta_xi_base_tempt;
    delete[] cross_xi_b;
    delete[] delta_eta_base_tempt;
    delete[] cross_xi_end;
    delete[] delta_xi_end_tempt;
    delete[] delta_eta_end_tempt;
        
}
