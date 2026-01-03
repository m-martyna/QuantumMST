import numpy as np
from qiskit.quantum_info import SparsePauliOp

#https://quantum.cloud.ibm.com/learning/en/courses/quantum-diagonalization-algorithms/sqd-overview

matrix = np.loadtxt('n_scatterer_1D/results/secular_matrix.txt').view(complex)
print(matrix)

qubit_operator = SparsePauliOp.from_operator(matrix)
print(qubit_operator)