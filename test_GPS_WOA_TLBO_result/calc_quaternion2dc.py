import numpy as np


def calc_quaternion2dc(q0, q1, q2, q3):
    """
    四元素转方向余弦矩阵
    :param q0: 标量
    :param q1:
    :param q2:
    :param q3:
    :return:
    """

    A_b = np.zeros((3, 3))

    A_b[0, 0] = q0 ** 2 + q1 ** 2 - q2 ** 2 - q3 ** 2
    A_b[0, 1] = 2 * (q1 * q2 - q0 * q3)
    A_b[0, 2] = 2 * (q1 * q3 + q0 * q2)
    A_b[1, 0] = 2 * (q1 * q2 + q0 * q3)
    A_b[1, 1] = q0 ** 2 - q1 ** 2 + q2 ** 2 - q3 ** 2
    A_b[1, 2] = 2 * (q2 * q3 - q0 * q1)
    A_b[2, 0] = 2 * (q1 * q3 - q0 * q2)
    A_b[2, 1] = 2 * (q2 * q3 + q0 * q1)
    A_b[2, 2] = q0 ** 2 - q1 ** 2 - q2 ** 2 + q3 ** 2

    return A_b
