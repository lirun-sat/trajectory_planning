import numpy as np
from cross import cross


def calc_v(v_b, w_b, r, r_b, A_links_transform, N, p, q_dot):
    v = np.zeros((N, 3, 1))
    Ez = np.array([[0],
                   [0],
                   [1]])

    for i in range(N):
        v[i, :, :] = v_b + cross(w_b) @ (r[i, :, :] - r_b)
        for j in range(i + 1):
            v[i, :, :] += cross(A_links_transform[j, :, :] @ Ez) @ (r[i, :, :] - p[j, :, :]) * q_dot[j, :]

    return v



