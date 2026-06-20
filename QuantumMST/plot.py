import numpy as np
import scipy.integrate
import string

import matplotlib.pyplot as plt
from .wavefunction import calc_psi
from .protein import protein_structure

plt.rcParams["mathtext.fontset"] = "stix"

def plot_potential(array_x, V, structure: protein_structure, title = None, path=None, title_size = 13, axis_size = 11):

        title_font = {
                'fontname': 'STIXGeneral',
                'size': title_size,
                'weight': 'semibold',
                'color': 'black'
        }

        axis_font = {
                'fontname': 'STIXGeneral',
                'size': axis_size,    
        }

        letters = string.ascii_uppercase

        V_max = np.max(V)
        V_min = np.min(V)
        dV = 2
        V = np.array(V)
        w = len(structure.scatterers)
        if(w>10): plt.figure(figsize=(w*2, 10))
        plt.plot(array_x, V, color = "black")
        # plt.fill_between(array_x, V, V_max, where=V<0, color='red', alpha=0.1, hatch="//")
        plt.title(title, fontdict=title_font)
        plt.xlabel("x[$a_0$]", fontdict=axis_font)
        plt.ylabel("V[Ry]", fontdict=axis_font)
        plt.ylim(V_min-dV, V_max+dV)
        plt.xlim(min(array_x), max(array_x))

        is_neg = V < 0
        changes = [i for i in range(1, len(V)) if is_neg[i] != is_neg[i - 1]]
        intervals = [0] + changes + [len(V)]
        scatterer_count = 0
        letters = string.ascii_uppercase

        for i in range(len(intervals) - 1):
                start, end = intervals[i], intervals[i + 1]
                x_mid = (array_x[start] + array_x[end - 1]) / 2
                
                if V[(start + end) // 2] < 0:
                        label = f"I{letters[scatterer_count]}"
                        scatterer_count += 1
                        plt.text(x_mid, -0.5, label, color="red", fontsize=12, ha="center", va="center", weight="bold")
                        plt.fill_between(array_x[start:end], V[start:end], 0,color='red', alpha=0.1, hatch="//")
                else:
                        plt.text(x_mid, -0.5, "II", color="black", fontsize=10, ha="center", va="center", weight="semibold")



        
        if(path is not None): plt.savefig(path)
        plt.show()
        plt.close()


def plot_wavefunction(array_x, V, structure: protein_structure):
        title_font = {
        'fontname': 'STIXGeneral',
        'size': 13,
        'weight': 'semibold',
        'color': 'black'
        }

        axis_font = {
                'fontname': 'STIXGeneral',
                'size': 13,    
        }
        psi = np.array([calc_psi(x, structure) for x in array_x])

        integral = scipy.integrate.simpson((np.abs(psi)**2), x = array_x)
        norm_psi = psi/np.sqrt(integral)

        a = np.max(np.abs(norm_psi))*1.1
        # print(f"a = {a}")
        V = np.array(V)
        plt.plot(array_x, norm_psi, color='black')
        plt.axhline(0, ls='--', color='black')
        plt.fill_between(array_x, -a, a, where=V<0, color='red', alpha=0.1, hatch="//")
        plt.ylim(-a, a)
        plt.title(f"Wavefunction for energy {structure.energy:4f} Ry \nfor two-scatterer system", fontdict=title_font)
        plt.xlabel("x[$a_0$]", fontdict=axis_font)
        plt.ylabel(r"$\Psi$", fontdict=axis_font)
        plt.show()