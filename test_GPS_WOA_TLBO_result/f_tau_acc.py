from calc_binomial import calc_binomial


def f_tau_acc(tau):
    return 5 * calc_binomial(4, 2) * (-2 * (1 - tau) * tau ** 2 + (1 - tau) ** 2 * 2 * tau)
