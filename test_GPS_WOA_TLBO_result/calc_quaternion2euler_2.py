import numpy as np
from math import atan2
from math import asin
from math import pi


def calc_quaternion2euler_2(q0, q1, q2, q3):
    """
    function: 四元素转欧拉角
    :param q0:
    :param q1:
    :param q2:
    :param q3:
    :return: RPY 角： 绕固定坐标系的 X-Y-Z 依次旋转  α, β, γ  角
    """

    Epsilon = 0.0009765625
    Threshold = 0.5 - Epsilon
    TEST = q0 * q2 - q1 * q3
    if TEST < -Threshold or TEST > Threshold:
        sign_of_TEST = np.sign(TEST)
        gamma = (-2) * sign_of_TEST * atan2(q1, q0)
        beta = sign_of_TEST * (pi / 2.0)
        alfa = 0
    else:
        alfa = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 ** 2 + q2 ** 2))
        beta = asin(2 * (q0 * q2 - q1 * q3))
        gamma = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 ** 2 + q3 ** 2))

    return alfa, beta, gamma




