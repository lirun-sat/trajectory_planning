import numpy as np


def cross(u):
    n = np.zeros((3, 3))
    n[0, 1] = -u[2, :]
    n[0, 2] = u[1, :]
    n[1, 2] = -u[0, :]
    n[1, 0] = -n[0, 1]
    n[2, 0] = -n[0, 2]
    n[2, 1] = -n[1, 2]

    return n
