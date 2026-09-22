"""
Core module for solving and evaluating the 1D quantum wave function (Psi) 
within a multiple-scattering formalism.

"""
from .protein import protein_structure
from .scatterer import scatterer
import numpy as np
import scipy.linalg as la
import scipy.integrate
import matplotlib.pyplot as plt

def j_1d(l, z):
    return np.cos(z) if l == 0 else np.sin(z)

def h_1d(l, z):
    return np.exp(1j * z) if l == 0 else -1j * np.exp(1j * z)

def calc_R_k(E, V, r, l):
    k = np.sqrt(E - V + 0j)
    z = k * r
    R = np.cos(z) if l == 0 else np.sin(z)
    return R, k

def calculate_b(matrix):
    A_hermitian = matrix.conj().T
    Q, R_qr, P = la.qr(A_hermitian, pivoting=True)
    
    tol = 10e-1
    rank = np.sum(np.abs(np.diag(R_qr)) > tol)
    dependent_indices = sorted(P[rank:])

    if len(dependent_indices) > 0:
        index = dependent_indices[0]
    else:
        index = 0

    matrix_1 = np.delete(matrix, index, axis=0)
    b_vec = -matrix_1[:, index]
    matrix_remaining = np.delete(matrix_1, index, axis=1)

    try:
        remaining_solutions = np.linalg.solve(matrix_remaining, b_vec)
        full_solution = np.insert(remaining_solutions, index, 1.0)
        return full_solution
    except np.linalg.LinAlgError:
        return np.ones(b_vec.size + 1)

def Psi_components(structure):
    sec_matrix = structure.secular_matrix
    scatterer_num = len(structure.scatterers)
    energy = structure.energy


    B = calculate_b(sec_matrix)
    C = np.zeros_like(B, dtype=complex)

    index = 0
    for i_sc in range(scatterer_num):
        s = structure.scatterers[i_sc]
        V_sc = s.V
        r_sc = s.r
        

        k_ext = np.sqrt(energy + 0j)
        z_ext = k_ext * r_sc
        
        for l in [0, 1]:
            R_loc, _ = calc_R_k(energy, V_sc, r_sc, l)
            

            m_l = s.calc_ml(energy, l) if hasattr(s, 'calc_ml') else s.m[l]

            numerator = B[index] * (h_1d(l, z_ext) + m_l * j_1d(l, z_ext))
            C[index] = 0.0 + 0j if np.abs(R_loc) < 1e-12 else numerator / R_loc
            index += 1

    return B, C

def region(x, structure):
    scatterer_num = len(structure.scatterers)
    for i in range(scatterer_num):
        x_s = structure.positions[i]
        r = structure.scatterers[i].r
        if np.abs(x - x_s) <= r:
            return i
    return -1

def Psi(x, structure):
    B, C = Psi_components(structure)
    scatterer_num = len(structure.scatterers)
    energy = structure.energy
    k_ext = np.sqrt(energy + 0j)
    
    reg = region(x, structure)
    
    if reg < 0:
        psi_val = 0.0 + 0j
        index = 0
        for i_sc in range(scatterer_num):
            dx = x - structure.positions[i_sc]
            r_n = np.abs(dx)
            z = k_ext * r_n
            
            for l in [0, 1]:
                Y = 1.0 / np.sqrt(2) if l == 0 else (np.sign(dx) if dx != 0 else 1.0) / np.sqrt(2)
                psi_val += B[index] * h_1d(l, z) * Y
                index += 1
        return psi_val
    else:
        psi_val = 0.0 + 0j
        dx = x - structure.positions[reg]
        r_n = np.abs(dx)
        s = structure.scatterers[reg]
        
        for l in [0, 1]:
            index = reg * 2 + l
            R_loc, _ = calc_R_k(energy, s.V, r_n, l)
            Y = 1.0 / np.sqrt(2) if l == 0 else (np.sign(dx) if dx != 0 else 1.0) / np.sqrt(2)
            psi_val += C[index] * R_loc * Y
            
        return psi_val
