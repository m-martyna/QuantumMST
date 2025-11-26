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

def random_x(n, a, b, min_dist):
    x = []
    while len(x) < n:
        p = random.uniform(a, b)
        if all(abs(p - i) >= min_dist for i in x):
            x.append(p)

    return np.sort(x)

def generate_scatterer(number, min_dist_pot, rmin, rmax,  Vmin, Vmax, xmin, xmax):
    s = []
    x = random_x(number, xmin, xmax, 2*rmin+min_dist_pot)
    r = 0

    for i in range(number):
        if(i==0): beg = xmin
        else: beg = x[i-1] + r + min_dist_pot

        if(i==len(x)-1): end = xmax
        else: end = x[i+1] - min_dist_pot - rmin

        dist = min(np.abs(x[i]-beg), np.abs(x[i]-end))

        r = random.uniform(rmin, min(rmax, dist))
        V = random.uniform(Vmin, Vmax)
        s.append(scatterer(r, V))

    return s, x

        