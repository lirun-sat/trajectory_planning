from rpy2dc import rpy2dc
import numpy as np


def calc_links_transform(A_b, rpy_joints, q, N):
    """
    计算每个连杆坐标系相对于惯性系的转换矩阵
    :param A_b: 基座姿态，用方向余弦矩阵表示
    :param rpy_joints: 各个关节坐标系相对关系，用RPY表示
    :param q: 关节角度变量
    :param N:
    :return:
    """
    A_links_transform = np.zeros((N, 3, 3))
    for i in range(N):
        if i == 0:
            tempt = A_b @ rpy2dc(rpy_joints[0, i], rpy_joints[1, i], rpy_joints[2, i])
            A_links_transform[i, :, :] = tempt @ rpy2dc(0, 0, q[i, :])
        else:
            A_links_transform[i, :, :] = A_links_transform[i - 1, :, :] @ rpy2dc(rpy_joints[0, i], rpy_joints[1, i], rpy_joints[2, i])
            A_links_transform[i, :, :] = A_links_transform[i, :, :] @ rpy2dc(0, 0, q[i, :])

    return A_links_transform

