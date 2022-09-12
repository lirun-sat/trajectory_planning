import numpy as np


def calc_r(r_b, b_b, A_b, A_links_transform, a, b, N):
    r = np.zeros((N, 3, 1))
    for i in range(N):
        if i == 0:
            r[i, :, :] = r_b + A_b @ b_b + A_links_transform[i, :, :] @ a[i, :, :]
        else:
            r[i, :, :] = r[i - 1, :, :] + A_links_transform[i - 1, :, :] @ b[i - 1, :, :] + A_links_transform[i, :, :] @ a[i, :, :]

    return r
