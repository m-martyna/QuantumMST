"""
Utility module for mapping secular scattering matrices into quantum operator 
representations (Pauli strings) and evaluating expectation values via an Estimator.

"""
import numpy as np
from qiskit.quantum_info import Pauli, SparsePauliOp

def decimal_to_base_4(number):

    # convert a decimal number to base 4 and return the digits as a list

    digits = []
    while number:
        digits.append(number % 4)
        number //= 4
    return list(reversed(digits))

def convert_matrix(secular_matrix):

    # convert the secular matrix to a SparsePauliOp

    num_qubits = np.round(np.log2(secular_matrix.shape[0])).astype(int)
    paulis = ['I', 'X', 'Y', 'Z']
    operator_size = 2**num_qubits
    labels = []
    coeffs = []
    for i in range(4**num_qubits):
        i_base4 = decimal_to_base_4(i)
        while len(i_base4) < num_qubits: i_base4.insert(0, 0)
        label = ''
        for number in i_base4: label = label + paulis[number]

        pauli_matrix=Pauli(label)
        coefficient = np.real(np.trace(np.dot(np.matrix(secular_matrix).getH(), pauli_matrix))/operator_size)

        labels.append(label)
        coeffs.append(coefficient)
    spo = SparsePauliOp(labels, coeffs)
    return spo


def check_size(matrix):

    # check if the size of the matrix is a power of 2

    size = matrix.shape[0]
    if (size & (size - 1)) != 0: 
        next_pow2 = 1 << (size - 1).bit_length()
        pad_amount = next_pow2 - size
        return np.pad(matrix, ((0, pad_amount), (0, pad_amount)), mode='constant', constant_values=99)
    return matrix
    
def cost_func(params, ansatz, hamiltonian, estimator):
    pub = (ansatz, [hamiltonian], [params])
    result = estimator.run(pubs=[pub]).result()
    energy = result[0].data.evs[0]
    #print(f"Cost: {energy}")
    return energy

