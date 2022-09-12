import numpy as np


def calc_p(r, A_links_transform, a, N):
    p = np.zeros((N, 3, 1))

    for i in range(N):
        p[i, :, :] = r[i, :, :] + A_links_transform[i, :, :] @ (-a[i, :, :])

    return p
