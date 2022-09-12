import numpy as np
from rpy2dc import rpy2dc


def calc_A_b(alpha_b, beta_b, gamma_b):
    A_b = rpy2dc(alpha_b, beta_b, gamma_b)

    return A_b
