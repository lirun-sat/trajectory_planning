import numpy as np
from scipy import linalg
from cross import cross


def calc_J_bm(m_b, m, r, r_b,  I_b_body, I_links_body, A_b, A_links_transform, a, N):
    """
    calc_J_bm
    :param m_b:
    :param m:
    :param r:
    :param r_b:
    :param I_b_body: 在基座本体系中表示的基座惯量
    :param I_links_body: 在各个连杆本体系中表示的连杆惯量
    :param A_b:
    :param A_links_transform:
    :param a:
    :param N:
    :return:
    """

    Ez = np.array([[0],
                   [0],
                   [1]])

    # 计算惯性系中表示的基座惯量
    I_b = A_b @ I_b_body @ A_b.T

    # 计算系统质心位置
    r_g = m_b * r_b
    for i in range(N):
        r_g = r_g + m[i] * r[i, :, :]
    r_g = r_g / (m_b + sum(m))

    # 计算各个关节的位置
    p = np.zeros((N, 3, 1))
    for i in range(N):
        p[i, :, :] = r[i, :, :] + A_links_transform[i, :, :] @ (-a[i, :, :])

    # 计算惯性系中表示的各个连杆的惯量
    I_links = np.zeros((N, 3, 3))
    for i in range(N):
        I_links[i, :, :] = A_links_transform[i, :, :] @ I_links_body[i, :, :] @ A_links_transform[i, :, :].T

    # equation 2-53
    Hw = I_b
    for i in range(N):
        Hw += I_links[i, :, :] + m[i] * cross(r[i, :, :] - r_b).T @ cross(r[i, :, :] - r_b)

    # equation 2-49
    JR = np.zeros((N, 3, N))
    for i in range(N):
        for j in range(i + 1):
            JR[i, :, j] = (A_links_transform[j, :, :] @ Ez).flatten()

    # equation 2-38
    JT = np.zeros((N, 3, N))
    for i in range(N):
        for j in range(i + 1):
            JT[i, :, j] = (cross(A_links_transform[j, :, :] @ Ez) @ (r[i, :, :] - p[j, :, :])).flatten()

    # equation 2-54
    Hwq = np.zeros((3, N))
    for i in range(N):
        Hwq += I_links[i, :, :] @ JR[i, :, :] + m[i] * cross(r[i, :, :] - r_b) @ JT[i, :, :]

    # equation 2-60
    Hs = (m_b + sum(m)) * cross(r_g - r_b) @ cross(r_g - r_b) + Hw

    # equation 2-37
    JTw = np.zeros((3, N))
    for i in range(N):
        JTw += m[i] * JT[i, :, :]

    # equation 2-61
    Hq = Hwq - cross(r_g - r_b) @ JTw

    # equation 2-62 & equation 2-63
    J_bm_w = (-1) * linalg.inv(Hs) @ Hq
    J_bm_v = (-1) * (JTw / (m_b + sum(m)) + cross(r_g - r_b) @ linalg.inv(Hs) @ Hq)

    return J_bm_w, J_bm_v
