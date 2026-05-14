import numpy as np
from scipy import linalg
from scipy.sparse.linalg import eigs

def NumpySolver(matrix):
    return np.real(np.linalg.eigvals(matrix))

def ScipySolver(matrix):
    return np.real(linalg.eigvals(matrix))

def ScipySparseSolver(matrix, k=2):
    results, _ = np.real(eigs(matrix, k=k))
    return results

def VQESolver(matrix):
    pass

def VQDSolver(matrix):
    pass

def NewSolver(matrix):
    pass