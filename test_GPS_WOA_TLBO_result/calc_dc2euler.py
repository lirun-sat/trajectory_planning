import numpy as np
import math


def calc_dc2euler(R):
    """
    :param R: 旋转矩阵
    :return: 欧拉角：先绕固定系 X 旋转，再绕固定系 Y 旋转，再绕固定系 Z 旋转  对应的角度
    """
    if R[2, 0] != 1 and R[2, 0] != -1:
        theta_y_1 = -math.asin(R[2, 0])
        theta_y_2 = math.pi - theta_y_1

        theta_x_1 = math.atan2(R[2, 1] / math.cos(theta_y_1), R[2, 2] / math.cos(theta_y_1))
        theta_x_2 = math.atan2(R[2, 1] / math.cos(theta_y_2), R[2, 2] / math.cos(theta_y_2))

        theta_z_1 = math.atan2(R[1, 0] / math.cos(theta_y_1), R[0, 0] / math.cos(theta_y_1))
        theta_z_2 = math.atan2(R[1, 0] / math.cos(theta_y_2), R[0, 0] / math.cos(theta_y_2))

    elif R[2, 0] == -1:
        theta_z_1 = 0
        theta_z_2 = 0

        theta_y_1 = math.pi / 2
        theta_y_2 = math.pi / 2

        theta_x_1 = theta_z_1 + math.atan2(R[0, 1], R[0, 2])
        theta_x_2 = theta_z_1 + math.atan2(R[0, 1], R[0, 2])

    else:
        theta_z_1 = 0
        theta_z_2 = 0
        theta_y_1 = -math.pi / 2
        theta_y_2 = -math.pi / 2
        theta_x_1 = -theta_z_1 + math.atan2(-R[0, 1], -R[0, 2])
        theta_x_2 = -theta_z_1 + math.atan2(-R[0, 1], -R[0, 2])

    return theta_z_1 * 180 / math.pi, theta_y_1 * 180 / math.pi, theta_x_1 * 180 / math.pi, theta_z_2 * 180 / math.pi, theta_y_2 * 180 / math.pi, theta_x_2 * 180 / math.pi

