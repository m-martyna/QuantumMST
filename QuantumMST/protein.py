import numpy as np
import random
from .scatterer import scatterer
from .struct_const import G_AB

class protein_structure:
    def __init__(self, energy, scatterers = None, positions = None, reference_system = None):
        
        self.energy = energy

        if reference_system is not None:
            scatterers, positions = example_potentials(reference_system)
            self.scatterers = scatterers
            self.positions = positions
            self.struct_const = G_AB(energy, positions)
        elif scatterers is not None and positions is not None:
            self.scatterers = scatterers
            self.positions = positions
            self.struct_const = G_AB(energy, positions)
        else:
            raise ValueError("Either reference_system or both scatterers and positions must be provided.")

    def secular_matrix(self):
        n = len(self.scatterers)
        secular_matrix = -self.struct_const.matrix_g()
        E = self.struct_const.E

        for i in range(n):
            secular_matrix[i, i] = self.scatterers[i].calc_m(E)

        secular_matrix = secular_matrix.transpose(0, 2, 1, 3).reshape(2*n, 2*n)
        return secular_matrix
    
    def potential(self, resolution = 100):
        positions = self.positions
        scatterers = self.scatterers
        
        i= np.argsort(positions)
        i_max = i[-1]
        i_min = i[0]
        r_max = scatterers[i_max].r
        r_min = scatterers[i_min].r


        dx = 2.0 # additinal space
        x_beg = positions[i_min]-r_min-dx
        x_end = positions[i_max]+r_max+dx
        array_x = np.linspace(x_beg, x_end, len(scatterers)*resolution)

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
    
def example_potentials(number):

    if number == 1:
        rA, VA, xA = 0.8, -5, -1
        rB, VB, RAB = 0.6, -6, 2

        s1 = scatterer(r=rA, V=VA)
        s2 = scatterer(rB, VB)

        scatterers=[s1, s2]
        positions=[xA, xA+RAB]

    elif number == 2:
        half_valley = [-2.0, -4.0, -6.0, -8.0, -10.0, -10.0, -8.0, -6.0, -4.0, -2.0]
        all_V = half_valley + half_valley
        all_r = [0.15] * len(all_V)
        
        scatterers = []
        positions = []
        current_x = -9.0
        
        for r_val, V_val in zip(all_r, all_V):
            scatterers.append(scatterer(r=r_val, V=V_val))
            positions.append(current_x)
            current_x += 0.9

    elif number == 3:
        side_depths = [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0, -9.5]
        all_V = side_depths + [-10.0] + side_depths[::-1]
        all_r = [0.15] * len(all_V)
        
        scatterers = []
        positions = []
        current_x = -9.5
        
        for r_val, V_val in zip(all_r, all_V):
            scatterers.append(scatterer(r=r_val, V=V_val))
            positions.append(current_x)
            current_x += 0.9

    elif number == 4:
        one_quarter = [-3.0, -6.0, -9.0, -6.0, -3.0, -1.0]
        half_profile = one_quarter + one_quarter[::-1]
        all_V = half_profile + half_profile[::-1]
        all_r = [0.12] * len(all_V)
        
        scatterers = []
        positions = []
        current_x = -10.0
        
        for r_val, V_val in zip(all_r, all_V):
            scatterers.append(scatterer(r=r_val, V=V_val))
            positions.append(current_x)
            current_x += 0.8
    else:
        raise ValueError("Example number not recognized. Please choose a valid example number.")
    
    return scatterers, positions

def random_x(n, a, b, min_dist):
    x = []
    while len(x) < n:
        p = random.uniform(a, b)
        if all(abs(p - i) >= min_dist for i in x):
            x.append(p)

    return np.sort(x)

def generate_scatterer(number, min_dist_pot, rmin, rmax,  Vmin, Vmax, xmin, xmax, save_path=None):
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

    if(save_path is not None):
        r_values = [obj.r for obj in s] 
        V_values = [obj.V for obj in s]
        data = np.column_stack((x, r_values, V_values))
        np.savetxt(f"{save_path}structure.txt", data, delimiter='\t')
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
