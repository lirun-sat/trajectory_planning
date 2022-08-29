#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void quaternion2euler(double w, double x, double y, double z, double* euler)
{
    double Epsilon = 0.0009765625;
    double Threshold = 0.5 - Epsilon;
    
    double TEST = w * y - x * z;
    
    double gamma;
    double beta;
    double alpha;
    
    if(TEST < -Threshold || TEST > Threshold)
    {
        int sign;
        if(TEST < 0)
            sign = (-1);
        else
            sign = (1);

        gamma = -2 * sign * (double)atan2(x, w); // yaw

        beta = sign * (M_PI / 2.0); // pitch

        alpha = 0; // roll
    }
    else
    {
        alpha = atan2(2 * (y * z + w * x), w * w - x * x - y * y + z * z);
        beta = asin(-2 * (x * z - w * y));
        gamma = atan2(2 * (x * y + w * z), w * w + x * x - y * y - z * z);
    }
    
    euler[0] = alpha;
    euler[1] = beta;
    euler[2] = gamma; // 绕固定轴转
    
    
}



/*
int main()
{
    double* euler = new double[3];
    
    quaternion2euler(0.653, -0.271, 0.653, 0.271, euler);
    
    cout << euler[0] * 180 / PI << "  " << euler[1] * 180 / PI << "  " << euler[2] * 180 / PI << endl;
    
    delete[] euler;
    
}

*/


