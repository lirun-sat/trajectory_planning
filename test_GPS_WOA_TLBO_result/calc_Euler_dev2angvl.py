import numpy as np
from math import sin
from math import cos


def calc_Euler_dev2angvl(alpha_z, beta_y, gamma_x, alpha_dot, beta_dot, gamma_dot):

    N_Phi = np.array([
        [0, -sin(alpha_z), cos(alpha_z) * cos(beta_y)],
        [0, cos(alpha_z), sin(alpha_z) * cos(beta_y)],
        [1, 0, -sin(beta_y)]
    ])
    Euler_dev = np.array([
        [alpha_dot],
        [beta_dot],
        [gamma_dot]
    ])
    omega = N_Phi @ Euler_dev

    return omega
