from calc_binomial import calc_binomial


def f_tau_vel(tau):
    return 5 * calc_binomial(4, 2) * (1 - tau) ** 2 * tau ** 2
