#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


double calc_q_ddot(double tau)
{
    return 5 * (calc_binomial(4, 2) * (-2) * (1-tau) * pow(tau, 2) + calc_binomial(4, 2) * pow((1-tau), 2) * 2 * tau);
}