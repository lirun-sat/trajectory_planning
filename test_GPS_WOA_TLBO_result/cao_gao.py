import numpy as np
from random import choice
import time
import math
from scipy import linalg




# q_initial = np.array([[0],
#                       [40 * math.pi / 180],
#                       [-27 * math.pi / 180],
#                       [0],
#                       [-5 * math.pi / 180],
#                       [0]])
#
# a = q_initial.shape[0]
# print(a)
#
# b = list(q_initial.flatten())
# print(b)

# cc = np.array([
#     [1, 2, 4],
#     [0, 7, 5],
#     [9, 6, 8]
# ])
#
# dd = cc[2, :]
# print(type(dd))
# print(dd)
# print(type(cc[:, 1]))
# print((cc[:, 1]))
#
# ee =[]
# ee.append(dd)
# ddd = cc[1, :]
# ee.append(ddd)
# print(ee)
# print(ee[0])



# minx = [2, 3, 1, 5, 4]
# aa =   [1, 2, 9, 6, 9]
# maxx = [4, 6, 8, 9, 8]
#
# for i in range(len(aa)):
#     if minx[i] <= aa[i] <= maxx[i]:
#         pass
#
#     elif aa[i] < minx[i]:
#         aa[i] = minx[i]
#
#     elif aa[i] > maxx[i]:
#         aa[i] = maxx[i]
#
# print(aa)


#
# from Bezier_calc_forward_kinematics import Bezier_calc_forward_kinematics
#
# P_i3 = np.array([
#     [0.2447],
#     [0.4904],
#     [0.3491],
#     [0.0970],
#     [-0.0407],
#     [-0.25397]
# ])
# delta_xi_end_delta_Pe_end, delta_eta_end, delta_xi_base, delta_eta_base, T_f, delta_locus_e = Bezier_calc_forward_kinematics(P_i3)
#
# Pe_desired = np.array([[1.7185],
#                        [0.1041],
#                        [1.7335]])
#
# # 1.77222
# # 0.1313
# # 1.75509
#
# print(delta_xi_end_delta_Pe_end)
# print(delta_eta_end)
#
# print(delta_xi_base)
# print(delta_eta_base)
#
# print(T_f)
# print(delta_locus_e)

from calc_quaternion2euler import calc_quaternion2euler
from euler_ZYX2dc import euler_ZYX2dc
from calc_zyx2quaternion import calc_zyx2quaternion


# q1 = 0.347785924506779
# q2 = 0.5187695013380103
# q3 = 0.06015368355712106
# q0 = 0.7786556938409857
#
# alpha_temp, beta_temp, gamma_temp = calc_quaternion2euler(q0, q1, q2, q3)
# print(alpha_temp * 180 / math.pi, beta_temp * 180 / math.pi, gamma_temp * 180 / math.pi)


# alpha = 35.92 * math.pi / 180
# beta = 45.6647 * math.pi / 180
# gamma = 40.68 * math.pi / 180

# A_temp = euler_ZYX2dc(alpha, beta, gamma)
# print(A_temp)

# q0, q1, q2, q3 = calc_zyx2quaternion(45 * math.pi / 180, -50 * math.pi / 180, -30 * math.pi / 180 )
# print(q0, q1, q2, q3)

# q0, q1, q2, q3 = calc_zyx2quaternion(0.719, -0.277, 1.474)
# print(q0, q1, q2, q3)

# alpha_temp, beta_temp, gamma_temp = calc_quaternion2euler(q0, q1, q2, q3)
# print(alpha_temp * 180 / math.pi, beta_temp * 180 / math.pi, gamma_temp * 180 / math.pi)
#


# a = [[1, 2], [3, 5], [8, 4, 2], [6, 0, 8, 2]]
# print(type(a))
# print(type(a[0]))
# print(a[0])
# a[0].append(23)
# print(a[0])
# print(min(a[2]))
#
# q_ddot_buff = [[], [], [], []]
#
# q_ddot_buff[0].append(1)
# q_ddot_buff[0].append(7)
# q_ddot_buff[0].append(4)
#
# q_ddot_buff[1].append(3)
# q_ddot_buff[1].append(2)
# q_ddot_buff[1].append(7)
#
# q_ddot_buff[2].append(6)
# print(q_ddot_buff)
#
# print(q_ddot_buff[0][0])
# print(q_ddot_buff[0][1])
# print(q_ddot_buff[0][2])
# print(q_ddot_buff[1][2])
#
# a = q_ddot_buff
# print(a[0])
#
# print(len(a))


#
# import matplotlib.pyplot as plt
#
#
# def hua_qiu(x, y, z, r, dense):
#     """
#         圆心坐标 半径 稠密程度
#     """
#     t = np.linspace(0, np.pi * 2, dense)
#     s = np.linspace(0, np.pi, dense)
#     t, s = np.meshgrid(t, s)             # 生成稠密网格点
#     x = x + r * np.sin(s) * np.cos(t)    # 球面坐标公式
#     y = y + r * np.sin(s) * np.sin(t)
#     z = z + r * np.cos(s)
#     return x, y, z
#
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# xx, yy, zz = hua_qiu(x=0, y=0, z=0, r=1, dense=100)
# ax.plot_surface(xx, yy, zz, rstride=1, cstride=1, cmap='gray', alpha=0.5)  # cmap='rainbow',
# plt.show()


























