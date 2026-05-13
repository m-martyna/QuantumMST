# **Quantum computing in MST**

This project aims to explore multiple scattering theory and its potential applications within the field of quantum computing. This research was conducted as part of my Master’s degree titled **_Application of quantum computing to electronic structure calculation using Green’s function method_**, which provides more complex analysis of the topics covered in this repository.

## Requirements
The code was developed and tested using the following versions of:
* **Python:** `3.12.2`
* **Qiskit:** `2.4.1`

### Installation
Ensure you have Python installed ([python.org](https://www.python.org/)). The required quantum computing framework can be set up via `pip`:
```bash
pip install qiskit==2.4.1
```
## Overview
### Multiple scattering theory
The core of calculations is based on W. H. Butler’s article
> **_Validity and accuracy of multiple-scattering theory_**  
> https://doi.org/10.1103/PhysRevB.41.2684.   

The scattering problem for a 1D system with disjoint, symmetric potentials is equivalent to solving the secular matrix eigenvalue problem. If a given energy results in a zero determinant for the system, it is considered an allowed energy level. This matrix is constructed using a list of **`scatterer`** objects initialized with the dimensions of each well, an array of the center positions of scatterers and a chosen energy.  

To verify the implementation of the matrix, the analytical results from the article were compared with numerical results obtained using the NumPy library. The potential (a), the determinant  (b) as well as the eigenvalues (c) for various energy values are plotted below.

<table>
  <tr>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/readme/Butler/potential.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/readme/Butler/determinant.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/readme/Butler/eigen_values.png" width="300" height="300" /></td>
  </tr>
  <tr>
    <td align="center">(a)</td>
    <td align="center">(b)</td>
    <td align="center">(c)</td>
  </tr>
</table>

The wave functions plotted for the determined energies (c)(d) perfectly match those provided in the article (a)(b).

<table>
  <tr>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/readme/Butler/Psi1_Butler.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/readme/Butler/Psi2_Butler.png" width="300" height="300" /></td>
  </tr>
    <tr>
    <td align="center">(a)</td>
    <td align="center">(b)</td>
  </tr>
  <tr>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/readme/Butler/Psi1.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/readme/Butler/Psi2.png" width="300" height="300" /></td>
  </tr>

  <tr>
    <td align="center">(c)</td>
    <td align="center">(d)</td>
  </tr>

</table>


### Quantum computing
A size of secular matrix (2n x 2n) is determined by a number of scatteres (n) that creates a system. The computational time for complex systems scales poorly as the number of objects increases.
<table>
  <tr>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/readme/complexity.png" width="500" height="300" /></td>
  </tr>
<table>
To improve efficiency rather than computing determinants, it is possible to find several eigenvalues and use interpolation to identify where they cross the zero. This problem description immediately brings to mind the Variational Quantum Eigensolver (VQE), one of the quantum algorithms designed to solve challenging eigenvalue problems.   
There is a few issues with using it. Firstly, most quantum algorithms are designed for Hermitian matrices, whereas the matrix representing this problem is complex symmetric. Even if we meet this requirement, VQE can only find the ground state eigenvalue. To ensure we are not missing any allowed energies, we need to calculate all of them.
   
In this study two approaches are explored:  
* using VQD algorithm, which allows caluclating higher eigen values,
* creating custom algorithm that employs energy as one of paramters in the cost function.

## Contact & attribution
**Author:** Martyna Migdałek    
**Supervised by**:  
* **Tomasz Stopa, PhD**
* **Jakub Haberko, PhD**  

**Institution:** AGH University of Krakow, Faculty of Physics and Applied Computer Science  
**Email:** migdalekm@student.agh.edu.pl
