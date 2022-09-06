#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"



void calc_para_range(double* result_min, double* result_max)
{
    
    double M1_temp;
    double M2_temp;
    
    double M1;
    double M2;
    
    double* result_1;
    double* result_1_temp;
    result_1 = new double[N];
    result_1_temp = new double[N];
    
    double* result_2;
    double* result_2_temp;
    result_2 = new double[N];
    result_2_temp = new double[N];
    
    M1_temp = pow((1-0.0001), 5) + 5 * pow((1-0.0001), 4) * 0.0001 + 10 * pow((1-0.0001), 3) * pow(0.0001, 2);
    M2_temp = 10 * pow((1-0.0001), 2) * pow(0.0001, 3) + 5 * (1-0.0001) * pow(0.0001, 4) + pow(0.0001, 5);
             
    for(int i = 0; i < N; i++)
    {
        result_1[i] = (joint_angle_min_limit_rad[i] - q_INITIAL[i] * M1_temp) / M2_temp;
        result_2[i] = (joint_angle_max_limit_rad[i] - q_INITIAL[i] * M1_temp) / M2_temp;
        
    }
        
    for(double tau = 0.001; tau < 1.0; tau += 0.0001)
    {
        M1 = pow((1 - tau), 5) + 5 * pow((1 - tau), 4) * tau + 10 * pow((1 - tau), 3) * pow(tau, 2);
        
        M2 = 10 * pow((1 - tau), 2) * pow(tau, 3) + 5 * (1 - tau) * pow(tau, 4) + pow(tau, 5);
        
        for(int i = 0; i < N; i++)
        {
            result_1_temp[i] = (joint_angle_min_limit_rad[i] - q_INITIAL[i] * M1) / M2;
            result_2_temp[i] = (joint_angle_max_limit_rad[i] - q_INITIAL[i] * M1) / M2;
            
            if(result_1_temp[i] >= result_1[i])
                result_1[i] = result_1_temp[i];
                
            if(result_2_temp[i] <= result_2[i])
                result_2[i] = result_2_temp[i];
            
        } 
    
    }   
    
    for(int i = 0; i < N; i++)
    {
        result_min[i] = result_1[i];
        result_max[i] = result_2[i];
    }
    
    delete[] result_1;
    delete[] result_1_temp;
    delete[] result_2;
    delete[] result_2_temp;
    
}
        
      



