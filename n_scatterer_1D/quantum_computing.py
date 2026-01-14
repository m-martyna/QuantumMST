import numpy as np
import matplotlib.pyplot as plt
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import efficient_su2
from qiskit_algorithms import VQD
from qiskit_algorithms.optimizers import COBYLA
from qiskit_algorithms.state_fidelities import ComputeUncompute
from qiskit.primitives import StatevectorEstimator, StatevectorSampler # Dodano Sampler
from qiskit.quantum_info.analysis import Z2Symmetries

from scatterer import scatterer
from struc_const import G_AB
from functions import secular_matrix, interpolation

# optimizer = COBYLA()
# estimator = StatevectorEstimator()

# matrix = np.loadtxt('n_scatterer_1D/results/secular_matrix.txt').view(complex)

# matrix = (matrix + matrix.conj().T) / 2
# print(f'Matrix:\n{matrix}\n')
# print(f'Eigen values: \n{np.linalg.eigvals(matrix)}')

# qubit_operator = SparsePauliOp.from_operator(matrix)
# print(f'\n{qubit_operator}')

# print(f'\n Eigen values: {np.linalg.eigvals(matrix)}')

# op = SparsePauliOp.from_operator(-matrix)

# ansatz = efficient_su2(op.num_qubits, reps=3)

# vqe = VQE(estimator, ansatz, optimizer)
# result = vqe.compute_minimum_eigenvalue(op)
# print(-result.eigenvalue.real)

# rA = 0.8
# VA = -5
# xA = -1
# rB = 0.6
# VB = -6
# RAB = 2

# s1 = scatterer(r=rA, V=VA)
# s2 = scatterer(rB, VB)
# nE = 1

# s = [s1, s2]
# tab_x = [xA, xA+RAB]


# energies = np.linspace(-5.5, -2.5, 10)
# eigen_val = []
# for e in energies:
#     matrix = secular_matrix(s, G_AB(energy =  e, positions = tab_x))
#     matrix = (matrix + matrix.conj().T) / 2
#     qubit_operator = SparsePauliOp.from_operator(matrix)
#     op = SparsePauliOp.from_operator(-matrix)
#     symmetries = Z2Symmetries.find_z2_symmetries(op)
#     if not symmetries.is_empty():
#         tapered_ops = symmetries.taper(op)
#         reduced_op = tapered_ops[0]
#         op = reduced_op
#     print(op.num_qubits)
#     ansatz = efficient_su2(op.num_qubits, reps=3)
#     vqe = VQE(estimator, ansatz, optimizer)
#     result = vqe.compute_minimum_eigenvalue(op)
#     eigen_val.append(-result.eigenvalue.real)

# plt.scatter(energies, eigen_val)
# plt.show()



optimizer = COBYLA()
estimator = StatevectorEstimator()
sampler = StatevectorSampler() 
fidelity = ComputeUncompute(sampler) 

matrix = np.loadtxt('n_scatterer_1D/results/secular_matrix.txt').view(complex)
matrix = (matrix + matrix.conj().T) / 2
print(f'Matrix:\n{matrix}\n')
print(f'Eigen values: \n{np.linalg.eigvals(matrix)}')

op = SparsePauliOp.from_operator(-matrix)
ansatz = efficient_su2(op.num_qubits, reps=3)

vqd = VQD(estimator, fidelity, ansatz, optimizer, k=2)
result = vqd.compute_eigenvalues(op)
print(-result.eigenvalues.real)

rA = 0.8
VA = -5
xA = -1
rB = 0.6
VB = -6
RAB = 2

s1 = scatterer(r=rA, V=VA)
s2 = scatterer(rB, VB)
s = [s1, s2]
tab_x = [xA, xA+RAB]

nE = 20
energies = np.linspace(-5.5, -2.5, nE)
eigen_val = []
k = 2
for e in energies:
    matrix = secular_matrix(s, G_AB(energy=e, positions=tab_x))
    matrix = (matrix + matrix.conj().T) / 2
    
    op = SparsePauliOp.from_operator(-matrix)
    symmetries = Z2Symmetries.find_z2_symmetries(op)
    if not symmetries.is_empty():
        tapered_ops = symmetries.taper(op)
        reduced_op = tapered_ops[0]
        op = reduced_op
    
    ansatz = efficient_su2(op.num_qubits, reps=3)
    
    vqd = VQD(estimator, fidelity, ansatz, optimizer, k=k, betas=[3.0, 3.0])
    result = vqd.compute_eigenvalues(op)
    eigen_val.append(-result.eigenvalues.real)

eigen_val = np.array(eigen_val)
interpolation(nE, k, eigen_val, energies)