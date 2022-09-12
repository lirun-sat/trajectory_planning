import numpy as np
from cross import cross


def calc_J_bE(r_e, r_b):
    """
    equation 2-32
    :param r_e:
    :param r_b:
    :return:
    """
    J_bE = np.zeros((6, 6))
    J_bE[0:3, 0:3] = np.eye(3)
    J_bE[0:3, 3:6] = -cross(r_e - r_b)
    J_bE[3:6, 0:3] = np.zeros((3, 3))
    J_bE[3:6, 3:6] = np.eye(3)

    return J_bE
