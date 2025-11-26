import numpy as np
import matplotlib.pyplot as plt

from struc_const import G_AB
from functions import secular


def plot_potential(V, array_x):

    plt.plot(array_x, V, color = "black")
    plt.title("V(x)")
    plt.xlabel("x")
    plt.ylabel("V")

    V_max = np.max(V)
    V_min = np.min(V)
    dV = 2
    plt.ylim(V_min-dV, V_max+dV)

    plt.show()

def plot_eigen(s, tab_x, e1, e2, n, det = True, values = True):

    energies = np.linspace(e1, e2, n)
    det_t = []

    for E in energies:
        str_constans = G_AB(energy = E , positions = tab_x)
        det, eigen_values = secular(s, str_constans)
        det_t.append(det)
        if(values):plt.scatter(np.full(len(eigen_values), E), eigen_values, color = 'black', s = 2)
    if(values):plt.show()

    if(det):plt.plot(energies, det_t, color = 'green')
    plt.axhline(y=0, color='black', linestyle='--')
    plt.xlabel('energy')
    if(det):plt.show()