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

def secular(scatterers, struct_const, only_eigen_val = False):
    n = len(scatterers)
    secular_matrix = -struct_const.matrix_g()
    E = struct_const.E

    for i in range(n):
        secular_matrix[i, i] = scatterers[i].calc_m(E)

    secular_matrix = secular_matrix.transpose(0, 2, 1, 3).reshape(2*n, 2*n)
    if(only_eigen_val): return np.real(np.linalg.eigvals(secular_matrix))
    
    return np.real(np.linalg.det(secular_matrix)), np.real(np.linalg.eigvals(secular_matrix))

def secular_matrix(scatterers, struct_const, save = True):
    n = len(scatterers)
    secular_matrix = -struct_const.matrix_g()
    E = struct_const.E

    for i in range(n):
        secular_matrix[i, i] = scatterers[i].calc_m(E)

    secular_matrix = secular_matrix.transpose(0, 2, 1, 3).reshape(2*n, 2*n)
    if(save): np.savetxt('n_scatterer_1D/results/secular_matrix.txt', secular_matrix.view(float), delimiter='\t', fmt='%.6f')
    return secular_matrix

def eigen_val(s, tab_x, e1, e2, n, plot_det = True, plot_values = True, save = True):

    energies = np.linspace(e1, e2, n)
    det_t = []
    eigen_t = []
    with open('n_scatterer_1D/results/eigen_values.txt', 'w') as file:
        # start_time = time.perf_counter()
        for E in energies:
            str_constans = G_AB(energy = E , positions = tab_x)
            det, eigen_values = secular(s, str_constans)
            det_t.append(det)
            eigen_t.extend(np.sort(eigen_values))
            if(plot_values):plt.scatter(np.full(len(eigen_values), E), eigen_values, color = 'black', s = 2)
        # end_time = time.perf_counter()
        # print(f"time: {end_time-start_time:.6f} seconds")
    if(plot_values):
        plt.axhline(y=0, color='red', linestyle='--')
        plt.show()

    if(plot_det):
        plt.plot(energies, det_t, 'o-', color = 'green', markersize=2)
        plt.axhline(y=0, color='black', linestyle='--')
        plt.xlabel('energy')
        plt.show()
    if(save):
        np.savetxt('n_scatterer_1D/results/determinant.txt',
                    det_t, delimiter='\t')    
        np.savetxt('n_scatterer_1D/results/eigen_values.txt',
                   eigen_t, delimiter='\t')
        np.savetxt('n_scatterer_1D/results/energies.txt',
                    energies)
        
def random_x(n, a, b, min_dist):
    x = []
    while len(x) < n:
        p = random.uniform(a, b)
        if all(abs(p - i) >= min_dist for i in x):
            x.append(p)

    return np.sort(x)

def generate_scatterer(number, min_dist_pot, rmin, rmax,  Vmin, Vmax, xmin, xmax, save = True):
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

    if(save):
        r_values = [obj.r for obj in s] 
        V_values = [obj.V for obj in s]
        data = np.column_stack((x, r_values, V_values))
        np.savetxt('n_scatterer_1D/results/random_scatterer.txt', data, delimiter='\t')
    return s, x

def generate_even_scatterer(number):
    dis = 3
    x = np.arange(0, number)*dis
    s = [scatterer(0.5, -5) for _ in range(number)]
    return s, x


def load_scatterer(filename):
    data = np.loadtxt(filename)
    x, r, V = data.T
    s=[]
    for i in range(len(x)):
        s.append(scatterer(r[i], V[i]))
    
    return s, x

def interpolation(nE, values = None, energies = None):
    
    note = []

    wyniki = []
    if values is None:
        values = np.loadtxt('n_scatterer_1D/results/eigen_values.txt')

    n = int(len(values)/nE)
    values = np.array(values).reshape(nE, n)
    sign = []

    for i in range(n): 
        v_max = np.max(values[:, i])
        v_min = np.min(values[:, i])
        if(np.sign(v_max*v_min)==-1): sign.append(i)

    # print(sign)
    if energies is None:
        energies = np.loadtxt('n_scatterer_1D/results/energies.txt')


    for index in sign:
        plt.plot(energies, values[:, index], 'o-', markersize=2)

        sign_prev = values[0, index]
        for i in range(1, nE):
            sign_current = values[i, index]
            if(np.sign(sign_prev*sign_current)==-1):
                # print(f'{energies[i-1]:.4f} : {energies[i]:.4f}')
                start = i-1
                stop = i
                x = []
                y = []
                for j in range(i-2, i+2):
                    x.append(energies[j])
                    y.append(values[j, index])
                coeff = np.polyfit(x, y, 3)
                p = np.poly1d(coeff)
                result = p.roots
                start = min(energies[i-1], energies[i])
                end = max(energies[i-1], energies[i])
                result = result[(result >= start) & (result <= end)]
                wyniki.append(result.real)
                # print(result.real)
                break

        wiersz = {
            'Eigenvalue number': f"{sign[index]}",
            'In between the points': f'{energies[i-1]:.4f} : {energies[i]:.4f}', 
            'Interpolated value': result.real[0]      
        }
            
        note.append(wiersz)
    
    plt.axhline(y=0, color='black', linestyle='--')
    # plt.show()


    df = pd.DataFrame(note)
    # df = df.set_index('Eigenvalue number')
    df.index.name = None
    display(df)

    return wyniki







