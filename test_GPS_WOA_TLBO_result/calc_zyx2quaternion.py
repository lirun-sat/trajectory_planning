import numpy as np
from math import cos
from math import sin


def calc_zyx2quaternion(psi, theta, phi):
    """
    欧拉角转四元素，其中欧拉角用 ZYX 法表示
    :param psi: 表示 Z 轴转动量
    :param theta: 表示 Y 轴转动量
    :param phi: 表示 X 轴转动量
    :return: 四元素，q0标量，q123向量
    """

    q0 = cos(phi/2) * cos(theta/2) * cos(psi/2) + sin(phi/2) * sin(theta/2) * sin(psi/2)

    q1 = sin(phi/2) * cos(theta/2) * cos(psi/2) - cos(phi/2) * sin(theta/2) * sin(psi/2)

    q2 = cos(phi/2) * sin(theta/2) * cos(psi/2) + sin(phi/2) * cos(theta/2) * sin(psi/2)

    q3 = cos(phi/2) * cos(theta/2) * sin(psi/2) - sin(phi/2) * sin(theta/2) * cos(psi/2)

    return q0, q1, q2, q3
