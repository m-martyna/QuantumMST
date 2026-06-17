import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.sparse.linalg import eigs

from .protein import protein_structure
from .scatterer import scatterer
from .vqd import convert_matrix, check_size, decimal_to_base_4, efficient_su2

from qiskit_algorithms import VQD
from qiskit_ibm_runtime import Session, Estimator, Sampler
from qiskit_algorithms.state_fidelities import ComputeUncompute
from qiskit_algorithms.optimizers import COBYLA

def NumpySolver(matrix):
    return np.real(np.linalg.eigvals(matrix))

def ScipySolver(matrix):
    return np.real(linalg.eigvals(matrix))

def ScipySparseSolver(matrix, k=2):
    results, _ = np.real(eigs(matrix, k=k))
    return results

def Numpy_energies(energies, arr_s, arr_x):
    results = []
    for e in energies:
        structure = protein_structure(e, scatterers=arr_s, positions=arr_x)
        matrix = structure.secular_matrix()

        eigen_values = np.sort(NumpySolver(matrix))

        results.append(eigen_values)
    
    return np.array(results).T

# def VQDSolver(backend, k, energies, scatterers, positions, b = 1.5):
#     nE = len(energies)
#     last_optimal_points = None

#     evals = []
#     for _ in range(k):
#         evals.append([])

#     with Session(backend=backend) as session:
#         for energy in energies:
#             structure = protein_structure(energy, scatterers=scatterers, positions=positions)
#             matrix = structure.secular_matrix()

#             matrix_neg = -matrix
#             matrix_neg = check_size(matrix_neg)
#             op = convert_matrix(matrix_neg)

#             ansatz = efficient_su2(op.num_qubits, reps=4).decompose()
#             optimizer = COBYLA(maxiter = 1000)
            
#             estimator = Estimator(mode=session)
#             sampler = Sampler(mode=session)
#             fidelity = ComputeUncompute(sampler)
            
#             vqd = VQD(estimator, fidelity, ansatz, optimizer, k=k, betas=[b]*k, initial_point=last_optimal_points)

#             result_vqd = vqd.compute_eigenvalues(op)

#             last_optimal_points = result_vqd.optimal_points

#             for ki in range(k):
#                 evals[ki].append(-result_vqd.eigenvalues[ki].real)
#     return evals

def NewSolver(matrix):
    pass
