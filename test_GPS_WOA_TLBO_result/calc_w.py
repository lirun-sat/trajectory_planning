import numpy as np


def calc_w(w_b, A_links_transform, N, q_dot):
    w = np.zeros((N, 3, 1))
    Ez = np.array([[0],
                   [0],
                   [1]])
    for i in range(N):
        if i == 0:
            w[i, :, :] = w_b + A_links_transform[i, :, :] @ Ez * q_dot[i, :]
        else:
            w[i, :, :] = w[i - 1, :, :] + A_links_transform[i, :, :] @ Ez * q_dot[i, :]

    return w

