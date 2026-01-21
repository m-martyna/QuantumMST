import numpy as np
import matplotlib.pyplot as plt
from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit.library import efficient_su2
from qiskit_algorithms import VQD
from qiskit_algorithms.optimizers import COBYLA
from qiskit_algorithms.state_fidelities import ComputeUncompute
from qiskit.primitives import StatevectorEstimator, StatevectorSampler
from qiskit.quantum_info.analysis import Z2Symmetries

from scatterer import scatterer
from struc_const import G_AB
from functions import secular_matrix, interpolation

from qiskit.quantum_info import DensityMatrix
from qiskit import QuantumCircuit



def xor(matrix):
    N = matrix.shape[0]
    num_qubits = int(np.log2(N))
    
    rho = DensityMatrix(matrix)
    qc = QuantumCircuit(num_qubits)
    
    qc.x(0)
    
    return rho.evolve(qc).data

def fwht(matrix):
    N = matrix.shape[0]
    n = int(np.log2(N))
    
    rho = DensityMatrix(matrix)
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    
    return rho.evolve(qc).data

def Pauli(matrix):
    N = matrix.shape[0]
    n = int(np.log2(N))
    A = matrix.copy()
    pauli_sparse = SparsePauliOp.from_operator(A)
    return pauli_sparse.simplify()

optimizer = COBYLA()
estimator = StatevectorEstimator()
sampler = StatevectorSampler() 
fidelity = ComputeUncompute(sampler)

s1 = scatterer(0.5, -3)
s2 = scatterer(0.8, -4)

energy = -3
s = [s1, s1]
tab_x = [-2, 0]
matrix1 = secular_matrix(s, G_AB(energy, positions=tab_x), save = False)
matrix1_xor = xor(matrix1)
matrix1_hadamard = fwht(matrix1_xor)
pauli1_result = Pauli(matrix1_hadamard)

op = pauli1_result
op = (op + op.adjoint()) / 2
ansatz = efficient_su2(op.num_qubits, reps=3)

vqd = VQD(estimator, fidelity, ansatz, optimizer, k=4)
result = vqd.compute_eigenvalues(op)
print(result.eigenvalues.real)
