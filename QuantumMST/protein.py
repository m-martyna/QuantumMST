"""
Manages multi-scatterer systems (protein structures),
calculates the secular matrix and potential profiles,
and provides synthetic/periodic potential generators.
"""


import numpy as np
import os
import random
from .scatterer import scatterer
from .struct_const import G_AB


class protein_structure:
    def __init__(self, energy = None, scatterers = None, positions = None, reference_system = None):
        if reference_system is not None:
            scatterers, positions = example_potentials(reference_system)
            self.scatterers = scatterers
            self.positions = positions

        elif scatterers is not None and positions is not None:
            self.scatterers = scatterers
            self.positions = positions

        else:
            raise ValueError("Either reference_system or both scatterers and positions must be provided.")

        self._energy = None
        self.struct_const = None
        self._secular_matrix = None
        # self.V = None
        # self.array_x = None
        self.potential(resolution=100, dx=2.0)

        if energy is not None:
            self.energy = energy

    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self, value):
        self._energy = value
        self.struct_const = G_AB(value, self.positions)
        self.calc_secular_matrix()

    def calc_secular_matrix(self):
        if self.struct_const is None:
            raise ValueError("Energy is not defined. Please set .energy first.")

        n = len(self.scatterers)
        mat = -self.struct_const.matrix_g()
        E = self.struct_const.E

        for i in range(n):
            mat[i, i] = self.scatterers[i].calc_m(E)

        mat = mat.transpose(0, 2, 1, 3).reshape(2*n, 2*n)
        self._secular_matrix = mat
        return mat

    @property
    def secular_matrix(self):
        if self._secular_matrix is None:
            raise ValueError("Secular matrix has not been computed yet. Set .energy first.")
        return self._secular_matrix

    def potential(self, resolution = 100, dx = 2.0):
        positions = np.array(self.positions)
        scatterers = self.scatterers

        i = np.argsort(positions)
        i_max = i[-1]
        i_min = i[0]
        r_max = scatterers[i_max].r
        r_min = scatterers[i_min].r

        x_beg = positions[i_min] - r_min - dx
        x_end = positions[i_max] + r_max + dx
        array_x = np.linspace(x_beg, x_end, len(scatterers) * resolution)

        V = []
        index = 0
        n_scatterers = len(i)

        for x in array_x:
            if index >= n_scatterers:
                V.append(0.0)
                continue

            current = i[index]
            a = positions[current] - scatterers[current].r
            b = positions[current] + scatterers[current].r

            if (x >= a) and (x <= b):
                V.append(scatterers[current].V)
            else:
                V.append(0.0)

            if x > b:
                index += 1

        self.V = np.array(V)
        self.array_x = array_x

        return self.V, self.array_x
    
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
    elif number == 5:
        all_V = [-1.2, -1.8, -3.5, -3.8, -5.5, -5.2, -7.4, -7.8, -9.3, -8.9, 
                 -10.0, 
                 -9.4, -8.7, -8.1, -6.8, -5.9, -5.1, -3.6, -2.8, -2.1, -0.9]
        
        all_r = [0.16, 0.25, 0.15, 0.19, 0.15, 0.22, 0.16, 0.17, 0.24, 0.15, 
                 0.23, 
                 0.15, 0.20, 0.16, 0.18, 0.15, 0.26, 0.16, 0.21, 0.16, 0.15]
        
        all_steps = [1.25, 0.85, 1.40, 0.95, 1.10, 0.75, 1.30, 0.90, 1.05, 1.50,
                     0.80, 1.20, 0.70, 1.35, 0.95, 1.15, 0.85, 1.45, 1.00, 1.10]
        
        scatterers = []
        positions = []
        current_x = -11.5
        
        for i, (r_val, V_val) in enumerate(zip(all_r, all_V)):
            scatterers.append(scatterer(r=r_val, V=V_val))
            positions.append(current_x)
            
            if i < len(all_steps):
                current_x += all_steps[i]

    elif number == 6:
        all_V = [-5.0] * 10
        all_r = [1.5] * 10
        all_steps = [4.5] * 9
        
        scatterers = []
        positions = []
        current_x = -20.25
        
        for i, (r_val, V_val) in enumerate(zip(all_r, all_V)):
            scatterers.append(scatterer(r=r_val, V=V_val))
            positions.append(current_x)
            
            if i < len(all_steps):
                current_x += all_steps[i]
                
    else:
        raise ValueError("Example number not recognized. Please choose a valid example number.")
    
    return scatterers, positions

def random_x(n, a, b, min_dist, max_attempts=10000):
    if (b - a) < (n - 1) * min_dist:
        raise ValueError(f"Interval [{a}, {b}] is too small to fit {n} points with a minimum separation of {min_dist}.")

    for _ in range(max_attempts):
        x = []
        attempts = 0
        while len(x) < n and attempts < 1000:
            p = random.uniform(a, b)
            if all(abs(p - item) >= min_dist for item in x):
                x.append(p)
            attempts += 1

        if len(x) == n:
            return np.sort(x)

    raise RuntimeError("Failed to sample points within the maximum number of attempts; constraints may be too strict or unsatisfiable.")

def generate_system(number, min_dist_pot, rmin, rmax, Vmin, Vmax, xmin, xmax, save_path=None):
    s = []
    safe_xmin = xmin + rmin
    safe_xmax = xmax - rmin
    min_dist_centers = 2 * rmin + min_dist_pot

    x = random_x(number, safe_xmin, safe_xmax, min_dist_centers)
    r = 0.0

    for i in range(number):
        if i == 0:
            beg = xmin
        else:
            beg = x[i - 1] + r + min_dist_pot

        if i == number - 1:
            end = xmax
        else:
            end = x[i + 1] - min_dist_pot - rmin

        max_allowed_r = min(rmax, x[i] - beg, end - x[i])
        max_allowed_r = max(rmin, max_allowed_r)

        r = random.uniform(rmin, max_allowed_r)
        V = random.uniform(Vmin, Vmax)
        s.append(scatterer(r, V))

    if save_path is not None:
        os.makedirs(save_path, exist_ok=True)
        r_values = [obj.r for obj in s]
        V_values = [obj.V for obj in s]
        data = np.column_stack((x, r_values, V_values))
        file_path = os.path.join(save_path, "structure.txt")
        np.savetxt(file_path, data, delimiter='\t', header="x\tr\tV", comments='')

    return s, x

def generate_periodic(number, dis=3.0, r=0.5, V=-5.0):
    if 2 * r >= dis:
        raise ValueError("Scatterer diameter (2*r) must be smaller than the lattice constant (dis).")
    
    x = np.arange(0, number) * dis
    s = [scatterer(r, V) for _ in range(number)]
    return s, x

def load_scatterer(filename):
    data = np.loadtxt(filename)
    x, r, V = data.T
    s=[]
    for i in range(len(x)):
        s.append(scatterer(r[i], V[i]))
    
    return s, x
