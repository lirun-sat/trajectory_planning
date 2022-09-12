import numpy as np
from scipy import linalg
import math
from calc_links_transform import calc_links_transform
from cross import cross
from calc_zyx2quaternion import calc_zyx2quaternion
from calc_quaternion2dc import calc_quaternion2dc
from calc_r import calc_r
from calc_p import calc_p
from calc_J_bm import calc_J_bm
from calc_J_bE import calc_J_bE
from f_tau_vel import f_tau_vel
from f_tau_acc import f_tau_acc


def Bezier_calc_forward_kinematics(P_i3):

    # <editor-fold desc=" system parameters ">
    N = 7
    Ez = np.array([[0],
                   [0],
                   [1]])

    rpy_joints = np.zeros((3, N))
    rpy_joints[:, 0] = np.array([0, 0, math.pi / 2])
    rpy_joints[:, 1] = np.array([-math.pi / 2, math.pi / 2, 0])
    rpy_joints[:, 2] = np.array([math.pi / 2, 0, math.pi])
    rpy_joints[:, 3] = np.array([-math.pi / 2, 0, 0])
    rpy_joints[:, 4] = np.array([math.pi / 2, 0, math.pi])
    rpy_joints[:, 5] = np.array([math.pi / 2, 0, 0])
    rpy_joints[:, 6] = np.array([-math.pi / 2, 0, 0])

    m_b = 1474.1
    b_b = np.array([[0],
                    [0],
                    [1.0864]])

    m = [39.618, 13.427, 22.757, 13.427, 42.614, 13.427, 8.269]

    a = np.zeros((N, 3, 1))
    a[0] = np.array([[0],
                     [0.013631],
                     [0.13304]])

    a[1] = np.array([[0],
                     [0.018673],
                     [0.08]])

    a[2] = np.array([[0],
                     [0.011017],
                     [0.8422]])

    a[3, :, :] = np.array([[0],
                           [0.018673],
                           [0.08]])

    a[4, :, :] = np.array([[0],
                           [0.11343],
                           [1.1109]])

    a[5, :, :] = np.array([[0],
                           [0.018673],
                           [0.08]])

    a[6, :, :] = np.array([[0],
                           [0],
                           [0.10565]])

    b = np.zeros((N, 3, 1))
    b[0, :, :] = np.array([[0],
                           [0.136369],
                           [0.056960]])

    b[1, :, :] = np.array([[0],
                           [0.1],
                           [0.08]]) \
                 - np.array([[0],
                             [0.018673],
                             [0.08]])

    b[2, :, :] = np.array([[0],
                           [0.1],
                           [1.08]]) \
                 - np.array([[0],
                             [0.011017],
                             [0.8422]])

    b[3, :, :] = np.array([[0],
                           [0.1],
                           [0.08]]) \
                 - np.array([[0],
                             [0.018673],
                             [0.08]])

    b[4, :, :] = np.array([[0],
                           [0.08],
                           [1.26]]) \
                 - np.array([[0],
                             [0.11343],
                             [1.1109]])

    b[5, :, :] = np.array([[0],
                           [0.1],
                           [0.08]]) \
                 - np.array([[0],
                             [0.018673],
                             [0.08]])

    b[6, :, :] = np.array([[0],
                           [0],
                           [0.16435]])
##########################################################################################################
    I_b_body = np.array([[17388.34, 0, 0],
                         [0, 1340.43, 0],
                         [0, 0, 17981.26]])

    I_links_body = np.zeros((N, 3, 3))
    I_links_body[0, :, :] = np.array([[0.30554, 0, 0],
                                      [0, 0.25403, 0.030757],
                                      [0, 0.030757, 0.15563]])

    I_links_body[1, :, :] = np.array([[0.043493, 0, 0],
                                      [0, 0.031896, 0],
                                      [0, 0, 0.029035]])

    I_links_body[2, :, :] = np.array([[2.69, 0, 0],
                                      [0, 2.68, 0],
                                      [0, 0, 0.06]])

    I_links_body[3, :, :] = np.array([[0.043493, 0, 0],
                                      [0, 0.031896, 0],
                                      [0, 0, 0.029035]])

    I_links_body[4, :, :] = np.array([[1.75, 0, 0],
                                      [0, 1.47, 0.29],
                                      [0, 0.29, 0.33]])

    I_links_body[5, :, :] = np.array([[0.043493, 0, 0],
                                      [0, 0.031896, 0],
                                      [0, 0, 0.029035]])

    I_links_body[6, :, :] = np.array([[0.04, 0, 0],
                                      [0, 0.04, 0],
                                      [0, 0, 0.01]])

    # </editor-fold>

    # q = np.array([[0],
    #               [0.78],
    #               [1.57],
    #               [0.78],
    #               [0],
    #               [-1.57],
    #               [0]])

    q = np.array([[0],
                  [0],
                  [0],
                  [0],
                  [0],
                  [0],
                  [0]])

    q_initial = q
    locus_e = 0  # locus 表示末端走过的路程

    # NEED TO BE MODIFIED
    joint_angle_velocity_min_limit = -0.0873
    joint_angle_velocity_max_limit = 0.0873
    joint_angle_acceleration_min_limit = -0.00873
    joint_angle_acceleration_max_limit = 0.00873

    rpy_base = np.array([0, 0, 0])
    Q_base_0, Q_base_1, Q_base_2, Q_base_3 = calc_zyx2quaternion(rpy_base[2], rpy_base[1], rpy_base[0])
    A_b = calc_quaternion2dc(Q_base_0, Q_base_1, Q_base_2, Q_base_3)
    eta_b = Q_base_0
    xi_b = np.array([[Q_base_1],
                     [Q_base_2],
                     [Q_base_3]])

    eta_b_initial = Q_base_0
    xi_b_initial = np.array([[Q_base_1],
                             [Q_base_2],
                             [Q_base_3]])

    A_links_transform = calc_links_transform(A_b, rpy_joints, q, N)

    M = m_b + sum(m)
    s1 = sum(m) * A_b @ b_b
    s2 = np.zeros((3, 1))
    for i in range(N):
        s2 = s2 + m[i] * A_links_transform[i, :, :] @ a[i, :, :]
    s3 = np.zeros((3, 1))
    for i in range(N):
        s3_temp = np.zeros((3, 1))
        for k in range(i):
            s3_temp = s3_temp + A_links_transform[k, :, :] @ (a[k, :, :] + b[k, :, :])
        s3 = s3 + m[i] * s3_temp
    r_b = -(s1 + s2 + s3) / M

    r = calc_r(r_b, b_b, A_b, A_links_transform, a, b, N)
    r_e = r[N-1, :, :] + A_links_transform[N-1, :, :] @ b[N-1, :, :]
    r_e_0 = r_e

    p = calc_p(r, A_links_transform, a, N)
    Jm_v = np.zeros((3, N))
    for i in range(N):
        Jm_v[:, i] = (cross(A_links_transform[i, :, :] @ Ez) @ (r_e - p[i, :, :])).flatten()

    Jm_w = np.zeros((3, N))
    for i in range(N):
        Jm_w[:, i] = (A_links_transform[i, :, :] @ Ez).flatten()

    Jm = np.concatenate((Jm_v, Jm_w), axis=0)
    J_bm_w, J_bm_v = calc_J_bm(m_b, m, r, r_b,  I_b_body, I_links_body, A_b, A_links_transform, a, N)
    J_bm = np.concatenate((J_bm_v, J_bm_w), axis=0)
    J_bE = calc_J_bE(r_e, r_b)
    J_g = J_bE @ J_bm + Jm
    J_g_v = J_g[0:3, :]
    J_g_w = J_g[3:6, :]

    # <editor-fold desc="  P_end_init and P_end_target  Q_end_init Q_end_target r_rot_direct r_rot_ang">
    Pe_init = r_e_0
    # print("Pe_init;", Pe_init)

    Pe_desired = np.array([[1.5],
                           [-1.4],
                           [2.2]])

    Pe_straight_line_locus = linalg.norm(Pe_desired - Pe_init)

    # rpy_end 表示末端相对于绝对坐标系的姿态，用RPY表示，这一步需要根据图2-5算一下
    rpy_end = np.array([0, 90 * math.pi / 180, 90 * math.pi / 180])
    Q_end_0, Q_end_1, Q_end_2, Q_end_3 = calc_zyx2quaternion(rpy_end[2], rpy_end[1], rpy_end[0])
    eta_end = Q_end_0
    xi_end = np.array([[Q_end_1], [Q_end_2], [Q_end_3]])

    Q_end_init = np.array([[Q_end_0], [Q_end_1], [Q_end_2], [Q_end_3]])

    ZYX_euler_end_desired_degree = np.array([45, -50, -30])
    ZYX_euler_end_desired_rad = ZYX_euler_end_desired_degree * math.pi / 180
    Q_end_0_desired, Q_end_1_desired, Q_end_2_desired, Q_end_3_desired = \
        calc_zyx2quaternion(ZYX_euler_end_desired_rad[0], ZYX_euler_end_desired_rad[1], ZYX_euler_end_desired_rad[2])

    eta_end_desired = Q_end_0_desired
    xi_end_desired = np.array([[Q_end_1_desired],
                               [Q_end_2_desired],
                               [Q_end_3_desired]])

    # </editor-fold>

    time_buff = []

    q_ddot_buff = [[], [], [], [], [], [], []]
    q_dot_buff = [[], [], [], [], [], [], []]
    q_buff = [[], [], [], [], [], [], []]

    r_e_buff = [[], [], []]
    Pe_buff = [[], [], []]
    eta_end_buff = []
    xi_end_buff = [[], [], []]
    eta_base_buff = []
    xi_base_buff = [[], [], []]
    Pe_initial_buff = []

    Pe = Pe_init
    Pe_initial_buff.append(Pe_init[0, 0])
    Pe_initial_buff.append(Pe_init[1, 0])
    Pe_initial_buff.append(Pe_init[2, 0])

    delta_xi_base_norm_max = 0

    # delta_tau = 0.001
    delta_tau = 0.01
    for tau in np.arange(0, 1, delta_tau):
        time_buff.append(tau)
        q_dot = f_tau_vel(tau) * (P_i3 - q_initial)
        q_ddot = f_tau_acc(tau) * (P_i3 - q_initial)
        q = q + q_dot * delta_tau
        for i in range(N):
            q_ddot_buff[i].append(q_ddot[i, :])

        for i in range(N):
            q_dot_buff[i].append(q_dot[i, :])

        for i in range(N):
            q_buff[i].append(q_dot[i, :])

        eta_b_dot = (- xi_b.T @ J_bm_w @ q_dot) / 2
        xi_b_dot = ((eta_b * np.eye(3) - cross(xi_b)) @ J_bm_w @ q_dot) / 2
        eta_end_dot = (- xi_end.T @ J_g_w @ q_dot) / 2
        xi_end_dot = ((eta_end * np.eye(3) - cross(xi_end)) @ J_g_w @ q_dot) / 2

        eta_b = eta_b + eta_b_dot * delta_tau
        xi_b = xi_b + xi_b_dot * delta_tau

        eta_b = eta_b / math.sqrt(eta_b * eta_b + xi_b[0, 0] * xi_b[0, 0] + xi_b[1, 0] * xi_b[1, 0] + xi_b[2, 0] * xi_b[2, 0])
        xi_b = xi_b / math.sqrt(eta_b * eta_b + xi_b[0, 0] * xi_b[0, 0] + xi_b[1, 0] * xi_b[1, 0] + xi_b[2, 0] * xi_b[2, 0])

        eta_end = eta_end + eta_end_dot * delta_tau
        xi_end = xi_end + xi_end_dot * delta_tau

        eta_end = eta_end / math.sqrt(eta_end * eta_end + xi_end[0, 0] * xi_end[0, 0] + xi_end[1, 0] * xi_end[1, 0] + xi_end[2, 0] * xi_end[2, 0])
        xi_end = xi_end / math.sqrt(eta_end * eta_end + xi_end[0, 0] * xi_end[0, 0] + xi_end[1, 0] * xi_end[1, 0] + xi_end[2, 0] * xi_end[2, 0])

        delta_xi_base = eta_b * xi_b_initial - eta_b_initial * xi_b - cross(xi_b) @ xi_b_initial
        delta_xi_base_norm = linalg.norm(delta_xi_base)
        if delta_xi_base_norm > delta_xi_base_norm_max:
            delta_xi_base_norm_max = delta_xi_base_norm

        # delta_eta_base = eta_b * eta_b_initial + xi_b.T @ xi_b_initial

        eta_end_buff.append(eta_end[0, 0])
        xi_end_buff[0].append(xi_end[0, 0])
        xi_end_buff[1].append(xi_end[1, 0])
        xi_end_buff[2].append(xi_end[2, 0])

        eta_base_buff.append(eta_b[0, 0])
        xi_base_buff[0].append(xi_b[0, 0])
        xi_base_buff[1].append(xi_b[1, 0])
        xi_base_buff[2].append(xi_b[2, 0])

        v_e = J_g_v @ q_dot
        Pe = Pe + v_e * delta_tau

        Pe_buff[0].append(Pe[0, 0])
        Pe_buff[1].append(Pe[1, 0])
        Pe_buff[2].append(Pe[2, 0])

        locus_e = locus_e + linalg.norm(v_e) * delta_tau

        v_b = J_bm_v @ q_dot
        r_b = r_b + v_b * delta_tau
        A_b = calc_quaternion2dc(eta_b, xi_b[0, :], xi_b[1, :], xi_b[2, :])
        A_links_transform = calc_links_transform(A_b, rpy_joints, q, N)
        r = calc_r(r_b, b_b, A_b, A_links_transform, a, b, N)
        r_e = r[N - 1, :, :] + A_links_transform[N - 1, :, :] @ b[N - 1, :, :]

        r_e_buff[0].append(r_e[0, 0])
        r_e_buff[1].append(r_e[1, 0])
        r_e_buff[2].append(r_e[2, 0])

        p = calc_p(r, A_links_transform, a, N)
        for i in range(N):
            Jm_v[:, i] = (cross(A_links_transform[i, :, :] @ Ez) @ (r_e - p[i, :, :])).flatten()
        for i in range(N):
            Jm_w[:, i] = (A_links_transform[i, :, :] @ Ez).flatten()

        Jm = np.concatenate((Jm_v, Jm_w), axis=0)
        J_bm_w, J_bm_v = calc_J_bm(m_b, m, r, r_b, I_b_body, I_links_body, A_b, A_links_transform, a, N)
        J_bm = np.concatenate((J_bm_v, J_bm_w), axis=0)
        J_bE = calc_J_bE(r_e, r_b)
        J_g = J_bE @ J_bm + Jm
        J_g_v = J_g[0:3, :]
        J_g_w = J_g[3:6, :]


    A_links_transform_EE = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            A_links_transform_EE[i, j] = A_links_transform[N - 1, i, j]
            A_links_transform_EE[i+3, j+3] = A_links_transform[N - 1, i, j]

    manipul = math.sqrt(linalg.det((A_links_transform_EE.T @ J_g) @ (A_links_transform_EE.T @ J_g).T))

    q_ddot_min_temp = []
    for i in range(N):
        q_ddot_min_temp.append(min(q_ddot_buff[i]))

    q_ddot_max_temp = []
    for i in range(N):
        q_ddot_max_temp.append(max(q_ddot_buff[i]))

    q_dot_min_temp = []
    for i in range(N):
        q_dot_min_temp.append(min(q_dot_buff[i]))

    q_dot_max_temp = []
    for i in range(N):
        q_dot_max_temp.append(max(q_dot_buff[i]))

    q_ddot_min = min(q_ddot_min_temp)
    q_ddot_max = max(q_ddot_max_temp)

    q_dot_min = min(q_dot_min_temp)
    q_dot_max = max(q_dot_max_temp)

    T_f = max(np.array([q_dot_min / joint_angle_velocity_min_limit,
                        q_dot_max / joint_angle_velocity_max_limit,
                        np.sqrt(q_ddot_min / joint_angle_acceleration_min_limit),
                        np.sqrt(q_ddot_max / joint_angle_acceleration_max_limit)]))

    # delta_xi_base = eta_b * xi_b_initial - eta_b_initial * xi_b - cross(xi_b) @ xi_b_initial
    # delta_eta_base = eta_b * eta_b_initial + xi_b.T @ xi_b_initial

    delta_xi_end = eta_end * xi_end_desired - eta_end_desired * xi_end - cross(xi_end) @ xi_end_desired
    delta_eta_end = eta_end * eta_end_desired + xi_end.T @ xi_end_desired
    delta_xi_end_norm = linalg.norm(delta_xi_end)

    delta_Pe_end = Pe_desired - Pe
    delta_Pe_end_norm = linalg.norm(delta_Pe_end)

    # delta_xi_end_delta_Pe_end = np.concatenate((delta_xi_end, delta_Pe_end), axis=0)

    # delta_xi_end_delta_Pe_end = delta_xi_end_delta_Pe_end.flatten()
    # delta_xi_base = delta_xi_base.flatten()
    #
    # delta_locus_e = locus_e - Pe_straight_line_locus

    return time_buff, \
           q_ddot_buff, q_dot_buff, q_buff, \
           Pe_buff, \
           r_e_buff, \
           eta_end_buff, xi_end_buff, \
           eta_base_buff, xi_base_buff, \
           Pe_initial_buff,  \
           r_e, eta_end, xi_end, \
           eta_b, xi_b, \
           manipul, T_f, locus_e, Pe_straight_line_locus, delta_xi_end_norm, delta_Pe_end_norm, delta_xi_base_norm_max
