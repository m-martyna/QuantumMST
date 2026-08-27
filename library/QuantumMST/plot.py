import numpy as np
import scipy.integrate
import string

import matplotlib.pyplot as plt
from .protein import protein_structure
from .wavefunction import Psi

plt.rcParams["mathtext.fontset"] = "stix"

def get_label_letters(n):
    result = []
    n += 1
    while n > 0:
        n, remainder = divmod(n - 1, 26)
        result.append(chr(65 + remainder))
    return "".join(reversed(result))

def plot_potential(structure: protein_structure, title = None, path=None, title_size = 13, axis_size = 11):

        V, array_x = structure.V, structure.array_x

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
                        label = f"I{get_label_letters(scatterer_count)}" #f"I{letters[scatterer_count]}"
                        scatterer_count += 1
                        plt.text(x_mid, -0.5, label, color="red", fontsize=12, ha="center", va="center", weight="bold")
                        plt.fill_between(array_x[start:end], V[start:end], 0,color='red', alpha=0.1, hatch="//")
                else:
                        plt.text(x_mid, -0.5, "II", color="black", fontsize=10, ha="center", va="center", weight="semibold")



        
        if(path is not None): plt.savefig(path)
        plt.show()
        plt.close()


def plot_wavefunction(array_x, V, structure: protein_structure):
    number_scatterers = len(structure.scatterers)
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

    psi = np.array([Psi(x, structure) for x in array_x])
    integral = scipy.integrate.simpson((np.abs(psi)**2), x = array_x)
    norm_psi = psi / np.sqrt(integral)

    a = np.max(np.abs(norm_psi)) * 1.1
    V = np.array(V)
    plt.plot(array_x, np.real(norm_psi), color='black')
    plt.axhline(0, ls='--', color='black')
    plt.fill_between(array_x, -a, a, where=V < 0, color='red', alpha=0.1, hatch="//")
    plt.ylim(-a, a)
    plt.title(f"Wavefunction for energy {structure.energy:.4f} Ry \nfor {number_scatterers}-scatterer system", fontdict=title_font)
    plt.xlabel("x[$a_0$]", fontdict=axis_font)
    plt.ylabel(r"$\Psi$", fontdict=axis_font)
    plt.show()