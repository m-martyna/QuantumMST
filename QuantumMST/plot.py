import numpy as np
import scipy.integrate 

import matplotlib.pyplot as plt
from .wavefunction import calc_psi
from .protein import protein_structure

def plot_potential(array_x, V, path=None):
        plt.plot(array_x, V, color = "black")
        plt.title("V(x)")
        plt.xlabel("x")
        plt.ylabel("V")
        V_max = np.max(V)
        V_min = np.min(V)
        dV = 2
        plt.ylim(V_min-dV, V_max+dV)
        if(path is not None): plt.savefig(path)
        plt.show()


def plot_wavefunction(array_x, V, structure: protein_structure):
        psi = np.array([calc_psi(x, structure) for x in array_x])

        integral = scipy.integrate.simpson((np.abs(psi)**2), x = array_x)
        norm_psi = psi/np.sqrt(integral)

        a = np.max(norm_psi)*1.1

        plt.plot(array_x, norm_psi, color='black')
        plt.axhline(0, ls='--', color='black')
        plt.fill_between(array_x, -a, a, where=V<0, color='red', alpha=0.1)
        plt.ylim(-a, a)
        # plt.title(f"E = {structure.energy}")
        plt.xlabel("x")
        plt.ylabel(f"$\Psi$")
        plt.show()