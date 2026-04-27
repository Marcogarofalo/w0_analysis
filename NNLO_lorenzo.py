import numpy as np
from scipy.special import kn
def g_tilda_1(lamb):
    j = np.array([1.,  2., 3., 4.,  5.,  6., 7.,  8.,  9.,], dtype=float)
    m = np.array([6., 12., 8., 6., 24., 24., 0., 12., 30.,], dtype=float)
    term = m * kn(1, lamb[:,np.newaxis] * np.sqrt(j))/np.sqrt(j)
    return 4*np.sum(term, axis=1)/lamb

def g_tilda_2(lamb):
    j = np.array([1.,  2., 3., 4.,  5.,  6., 7.,  8.,  9.,], dtype=float)
    m = np.array([6., 12., 8., 6., 24., 24., 0., 12., 30.,], dtype=float)
    term = m * kn(2, lamb[:,np.newaxis] * np.sqrt(j))/j
    return 4*np.sum(term, axis=1)/(lamb**2)


def Delta_pi_NNLO_CDH(L_a, xi, af_PS,f_or_M="f"):
    l_phy_1 = -0.4
    l_phy_2 = 4.3
    l_phy_3 = 3.2
    l_phy_4 = 4.4
    log_corr = 2*np.log(135.0/(np.sqrt(xi)*(af_PS*4*np.pi)))
    l_bar_1 = l_phy_1 + log_corr
    l_bar_2 = l_phy_2 + log_corr
    l_bar_3 = l_phy_3 + log_corr
    l_bar_4 = l_phy_4 + log_corr

    C_f_1 = -(7./9) + 2*l_bar_1 + (4./3)*l_bar_2 - 3*l_bar_4 + 4*l_bar_4 # last term is correction due to the use of measured xi instead of chi_PT prediction in the expansion
    C_f_2 = (112./9) - (8./3)*l_bar_1 - (32./3)*l_bar_2
    C_M_1 = -(55./18) + 4*l_bar_1 + (8./3)*l_bar_2 - (5./2)*l_bar_3 -2*l_bar_4 + 2*l_bar_4 # last term is correction due to the use of measured xi instead of chi_PT prediction in the expansion
    C_M_2 = C_f_2

    s_0 = 2-np.pi/2
    s_1 = np.pi/4 - 1./2
    s_2 = 1./2 - np.pi/8
    s_3 = 3*np.pi/16 - 1./2

    lamb = np.sqrt(xi) * 4 * np.pi * af_PS * L_a

    g1 = g_tilda_1(lamb)
    g2 = g_tilda_2(lamb)
    S_f_4 = (4./3*s_0 - 13./6*s_1)*g1 - (40/3*s_0 - 4*s_1 - 8./3*s_2 - 13./3*s_3)*g2
    S_M_4 = 13./3*s_0*g1 - (40./3*s_0 + 32./3*s_1 + 26./3*s_2)*g2

    if f_or_M == "f":
        return -2.0 * xi * g1 + 2 *xi**2 * (C_f_1 * g1 + C_f_2 * g2 + S_f_4)
    elif f_or_M == "M":
        return 0.5 * xi * g1 - xi**2 * (C_M_1 * g1 + C_M_2 * g2 + S_M_4)
    
x=Delta_pi_NNLO_CDH(1, np.array([1]), np.array([1]),f_or_M="f")
# Or using an f-string (cleaner)
print("FVE NNL ",f"{x[0]:.20f}")
x=Delta_pi_NNLO_CDH(1, np.array([1]), np.array([1]),f_or_M="M")
# Or using an f-string (cleaner)
print("FVE NNL ",f"{x[0]:.20f}")
