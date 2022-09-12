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
from f_tau_vel import f_tau_2
from f_tau_acc import f_tau_acc



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
#########################################################################################################
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

# q = np.array([[3.133],
#               [0.108],
#               [0.495],
#               [1.812],
#               [-2.595],
#               [-1.028],
#               [1.93]])

# q = np.array([[-0.484],
#               [0.227],
#               [0.819],
#               [2.030],
#               [0.262],
#               [1.043],
#               [2.111]])

q = np.array([[0],
                  [0],
                  [0],
                  [0],
                  [0],
                  [0],
                  [0]])

q_initial = q


rpy_base = np.array([0, 0, 0])
Q_base_0, Q_base_1, Q_base_2, Q_base_3 = calc_zyx2quaternion(rpy_base[2], rpy_base[1], rpy_base[0])
A_b = calc_quaternion2dc(Q_base_0, Q_base_1, Q_base_2, Q_base_3)

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



det_J_g_J_g_T = linalg.det(J_g @ J_g.T)

print(det_J_g_J_g_T)

