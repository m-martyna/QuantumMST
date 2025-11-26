import numpy as np
import random

from scatterer import scatterer

def potential(scatterers, positions):

    i= np.argsort(positions)
    i_max = i[-1]
    i_min = i[0]
    r_max = scatterers[i_max].r
    r_min = scatterers[i_min].r


    dx = 1.0
    x_beg = positions[i_min]-r_min-dx
    x_end = positions[i_max]+r_max+dx
    array_x = np.linspace(x_beg, x_end, len(scatterers)*100)

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


    return V, array_x

def secular(scatterers, struct_const):
    n = len(scatterers)
    secular_matrix = -struct_const.matrix_g()
    E = struct_const.E

    for i in range(n):
            secular_matrix[i, i] = scatterers[i].calc_m(E)

    secular_matrix = secular_matrix.transpose(0, 2, 1, 3).reshape(2*n, 2*n)

    return np.real(np.linalg.det(secular_matrix)), np.real(np.linalg.eigvals(secular_matrix))

def generate_scatteters(number, Vmin, Vmax, rmin, rmax):
    s = [scatterer(random.uniform(rmin, rmax), random.uniform(Vmin, Vmax)) for _ in range(number)]
    return s

def generate_positions(s, xmin, xmax):
    x = [random.uniform(xmin, xmax) for i in range(len(s))]
    return x