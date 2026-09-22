"""
Generates and plots
the potential profile V(x).
"""

import numpy as np
import random
import matplotlib.pyplot as plt
import time
import pandas as pd

from .scatterer import scatterer
from .struct_const import G_AB


def potential(scatterers, positions, plot = True):

    i= np.argsort(positions)
    i_max = i[-1]
    i_min = i[0]
    r_max = scatterers[i_max].r
    r_min = scatterers[i_min].r


    dx = 1.0
    x_beg = positions[i_min]-r_min-dx
    x_end = positions[i_max]+r_max+dx
    array_x = np.linspace(x_beg, x_end, len(scatterers)*500)

    V = []
    index = 0
    n_scatterers = len(i)

    for x in array_x:

        if index >= n_scatterers:
            V.append(0)
            continue
        
        current = i[index]
        a = positions[current]-scatterers[current].r
        b = positions[current]+scatterers[current].r

        if(x>=a) and (x<=b): V.append(scatterers[current].V)        
        else: V.append(0)

        if(x>b): index+=1

    if(plot):
        plt.plot(array_x, V, color = "black")
        plt.title("V(x)")
        plt.xlabel("x")
        plt.ylabel("V")
        V_max = np.max(V)
        V_min = np.min(V)
        dV = 2
        plt.ylim(V_min-dV, V_max+dV)
        plt.show()


    return V



