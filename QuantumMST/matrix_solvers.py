"""
Solvers and eigenvalue interpolation tools for finding bound states,
including classical numerical solvers (NumPy, FDM)
and quantum variational algorithms (VQD).
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import linalg
from scipy.sparse.linalg import eigs
from typing import Literal
from qiskit_aer import AerSimulator

from .protein import protein_structure
from .scatterer import scatterer
from .vqd import convert_matrix, check_size, decimal_to_base_4, efficient_su2

from qiskit_algorithms import VQD
from qiskit_ibm_runtime import Session, Estimator, Sampler
from qiskit_algorithms.state_fidelities import ComputeUncompute
from qiskit_algorithms.optimizers import COBYLA


title_font = {
        'fontname': 'STIXGeneral',
        'size': 13,
        'weight': 'semibold',
        'color': 'black'
}

axis_font = {
        'fontname': 'STIXGeneral',
        'size': 20,  
}


def interpolation(nE, values=None, energies=None, display_stats=False, plot=False):

    note = []

    wyniki = []
    if values is None:
        values = np.loadtxt("n_scatterer_1D/results/eigen_values.txt")

    n = int(len(values) / nE)
    values = np.array(values).reshape(nE, n)
    sign = []

    for i in range(n):
        v_max = np.max(values[:, i])
        v_min = np.min(values[:, i])
        if np.sign(v_max * v_min) == -1:
            sign.append(i)

    # print(sign)
    if energies is None:
        energies = np.loadtxt("n_scatterer_1D/results/energies.txt")

    if plot:
        plt.figure(figsize=(10, 6))
        cmap = plt.cm.coolwarm
        num_lines = len(sign)
        colors = [cmap(i) for i in np.linspace(0, 1, num_lines)]
        plt.ylabel('Eigenvalue',fontdict = axis_font)
        plt.xlabel('Energy[Ry]',fontdict = axis_font)
        norm = plt.Normalize(vmin=1, vmax=num_lines)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar1 = plt.colorbar(sm, ax=plt.gca(), ticks=[1, num_lines])
        cbar1.set_label(f'Eigenvalues count: {num_lines}', fontdict=axis_font)
        cbar1.ax.set_yticklabels(['Last', 'First'])
        plt.axhline(y=0, color="black", linestyle="--")

    c = 0
    for index in sign:
        if plot:
            plt.plot(energies, values[:, index], "o-", markersize=3, color=colors[c])
        c+=1
        sign_prev = values[0, index]
        for i in range(1, nE):
            sign_current = values[i, index]
            if np.sign(sign_prev * sign_current) == -1:
                # print(f'{energies[i-1]:.4f} : {energies[i]:.4f}')
                start0 = i - 1
                stop0 = i
                x = []
                y = []
                j_start = max(0, min(i - 2, nE - 4))
                j_end = min(nE, j_start + 4)
                for j in range(j_start, j_end):
                    x.append(energies[j])
                    y.append(values[j, index])
                coeff = np.polyfit(x, y, 3)
                p = np.poly1d(coeff)
                result = p.roots
                start = min(energies[i - 1], energies[i])
                end = max(energies[i - 1], energies[i])
                result = result[(result >= start) & (result <= end)]
                wyniki.append(result.real[0])
                # print(result.real)
                break
            sign_prev = sign_current

        wiersz = {
            "Eigenvalue number": f"{index}",
            "In between the points": f"{energies[start0]:.6f} : {energies[stop0]:.6f}",
            "Interpolated value": (
                result.real[0] if len(result) > 0 else np.nan
            ),
        }

        note.append(wiersz)

    
    if plot:
        plt.show()

    df = pd.DataFrame(note)
    # df = df.set_index('Eigenvalue number')
    df.index.name = None
    if display_stats:
        display(df)

    return wyniki

def Solvers(structure: protein_structure, mode , nE = 20, display_stats = False, plot = False, real = True, imaginary = False, backend = None, b = None, k = None):
    if mode == 'numpy':
        return NumpySolver(structure, nE, display_stats=display_stats, plot=plot, real=real, imaginary=imaginary)
    elif mode == 'vqd':
        return VQDSolver(structure, nE, backend, k, b, display_stats=display_stats, plot=plot)
    elif mode == 'new':
        return NewSolver(structure)
    else:
        raise ValueError("Invalid mode. Choose either 'numpy' or 'vqd'.")


def VQDSolver(
    structure, nE, backend=None, k=2, b=None, display_stats=False, plot=False):
    E_min = np.min(structure.V)
    energies = np.linspace(E_min, -0.1, nE)

    if backend is None:
        backend = AerSimulator()

    nE = len(energies)
    last_optimal_points = None
    arr_ev = []

    with Session(backend=backend) as session:
        for energy in energies:
            structure.energy = energy
            matrix = structure.secular_matrix

            matrix_neg = -matrix
            matrix_neg = check_size(matrix_neg)
            op = convert_matrix(matrix_neg)

            ansatz = efficient_su2(op.num_qubits, reps=4).decompose()
            optimizer = COBYLA(maxiter=1000)

            estimator = Estimator(mode=session)
            sampler = Sampler(mode=session)
            fidelity = ComputeUncompute(sampler)

            vqd = VQD(
                estimator=estimator,
                fidelity=fidelity,
                ansatz=ansatz,
                optimizer=optimizer,
                k=k,
                betas=b,
                initial_point=last_optimal_points,
            )

            result_vqd = vqd.compute_eigenvalues(op)
            last_optimal_points = result_vqd.optimal_points

            current_eigenvalues = -np.real(result_vqd.eigenvalues)
            current_eigenvalues = np.sort(current_eigenvalues)

            arr_ev.extend(current_eigenvalues)

    return np.array(interpolation(nE, arr_ev, energies, display_stats=display_stats, plot=plot))

def NumpySolver(structure: protein_structure, nE, display_stats = False, plot = False, real = True, imaginary = False):
    E_min = np.min(structure.V)
    energies = np.linspace(E_min, -0.1, nE) 
    arr_ev = []

    for E in energies:
        structure.energy = E
        matrix = structure.secular_matrix
        
        if real and imaginary:
            eigenvalue = np.sort(np.real(np.linalg.eigvals(matrix))) + 1j * np.sort(np.imag(np.linalg.eigvals(matrix)))
        elif imaginary:
            eigenvalue = np.sort(np.imag(np.linalg.eigvals(matrix)))
        else:
            eigenvalue = np.sort(np.real(np.linalg.eigvals(matrix)))
        arr_ev.extend(eigenvalue)

    return np.array(interpolation(nE, arr_ev, energies, display_stats=display_stats, plot=plot))

def NewSolver(structure: protein_structure):
    print("The algorithm has not been deployed yet while its cost function form is being finalized.")
    pass
       
def Numpy_energies(energies, arr_s, arr_x, real = True, imaginary = False):
    results = []
    for e in energies:
        structure = protein_structure(e, scatterers=arr_s, positions=arr_x)
        matrix = structure.secular_matrix

        eigen_values = np.sort(np.real(np.linalg.eigvals(matrix)))

        if real and imaginary:
            eigen_values = np.sort(np.real(np.linalg.eigvals(matrix))) + 1j * np.sort(np.imag(np.linalg.eigvals(matrix)))
        elif imaginary:
            eigen_values = np.sort(np.imag(np.linalg.eigvals(matrix)))
        else:
            eigen_values = np.sort(np.real(np.linalg.eigvals(matrix)))

        results.append(eigen_values)
    
    return np.array(results).T

def determinant(energies, arr_s, arr_x):
    results = []
    for e in energies:
        structure = protein_structure(e, scatterers=arr_s, positions=arr_x)
        matrix = structure.secular_matrix()
        results.append(np.real(np.linalg.det(matrix)))
        
    
    return results


def FDM(system: protein_structure, resolution):
    V, x = system.potential(resolution=resolution)

    N = len(V)
    dx = x[1] - x[0]

    diag = 2.0 / (dx**2) + V
    off_diag = -1.0 / (dx**2) * np.ones(N - 1)

    num_states = 12
    energies, wavefunctions = linalg.eigh_tridiagonal(
        diag, off_diag, select='i', select_range=(0, num_states - 1)
    )

    bound_mask = energies < 0
    bound_energies = energies[bound_mask]
    bound_wavefunctions = wavefunctions[:, bound_mask]

    for n, E in enumerate(bound_energies):
        print(f"E_{n} = {E:8.4f} Ry")


