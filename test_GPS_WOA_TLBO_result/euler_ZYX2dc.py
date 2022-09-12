import numpy as np
from math import cos
from math import sin


# Z-Y-X 欧拉角[alpha, beta, gamma]定义为：frame_j分别绕其 Z 轴、旋转后的 Y 轴、再旋转后的 X 轴旋转alpha, beta, gamma角后, 与frame_i重合

def euler_ZYX2dc(alpha, beta, gamma):
    A = np.zeros((3, 3))
    A[0, 0] = cos(alpha) * cos(beta)
    A[0, 1] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma)
    A[0, 2] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma)

    A[1, 0] = sin(alpha) * cos(beta)
    A[1, 1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma)
    A[1, 2] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma)

    A[2, 0] = -sin(beta)
    A[2, 1] = cos(beta) * sin(gamma)
    A[2, 2] = cos(beta) * cos(gamma)

    return A
