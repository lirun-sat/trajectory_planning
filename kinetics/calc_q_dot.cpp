#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"



double calc_q_dot(double tau)
{
    // 根据 Beseir 曲线 计算得到
    return (5 * calc_binomial(4, 2) * pow((1 - tau), 2) * pow(tau, 2));
}
