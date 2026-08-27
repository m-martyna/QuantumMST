# __Quantum computing in Multiple Scattering Theory (MST)__

This project aims to explore multiple scattering theory and its potential applications within the field of quantum computing. This research was conducted as part of my Master’s degree titled **_Application of quantum computing to electronic structure calculation using Green’s function method_**, which provides more complex analysis of the topics covered in this repository.

## Getting Started
### Requirements
The code was developed and tested using the following versions of:
* **Python:** `3.12.2`
* **Qiskit:** `2.4.1`

### Installation
This library can be installed directly via `pip`:
```bash
pip install qmst
```

## Overview

### Repository Structure

```text
├── QuantumMST/
│   ├── __init__.py          # Package initialization
│   ├── functions.py         # Helper mathematical and interpolation functions
│   ├── matrix_solvers.py    # Secular matrix solvers (MST and numerical methods)
│   ├── plot.py              # Visualization routines
│   ├── protein.py           # Protein structure and system definitions
│   ├── scatterer.py         # Definitions and building routines for single scatterer data
│   ├── struct_const.py      # Structural constants computation
│   ├── vqd.py               # Functions required to run the VQD (Variational Quantum Deflation) algorithm
│   └── wavefunction.py      # Wavefunction calculation and construction
│
├── 1_Introduction_MST.ipynb # Example usage: introduction to MST
├── 2_Solvers.ipynb          # Example usage: energy calculation
└── 3_Wavefunction.ipynb     # Example usage: wavefunction plotting
```


* **`QuantumMST/`** – Core Python package containing all classes, functions, and utilities required for calculations, quantum algorithms and visualization.
* **`1_Introduction_MST.ipynb`, `2_Solvers.ipynb`, `3_Wavefunction.ipynb`** – Jupyter notebooks demonstrating step-by-step example usage of the library.

### Multiple scattering theory
The core of calculations is based on W. H. Butler’s article
> **_Validity and accuracy of multiple-scattering theory_**  
> https://doi.org/10.1103/PhysRevB.41.2684.   


The scattering problem for a 1D system with disjoint, symmetric potentials is equivalent to solving the secular matrix eigenvalue problem. If a given energy results in a zero determinant for the system, it is considered an allowed energy level. This matrix is constructed using a list of **`scatterer`** objects initialized with the dimensions of each well, an array of the center positions of scatterers and a chosen energy.  

To verify the implementation of the matrix, the analytical results from the article were compared with numerical results obtained using the NumPy library. The potential (a), the determinant  (b) as well as the eigenvalues (c) for various energy values are plotted below.

#### Potential profile, determinants and eigenvalues
<table>
  <tr>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/.thesis_images/systems/2scatterers.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/.thesis_images/determinants/det_2s.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/.thesis_images/eigenvalues/ev2.png" width="300" height="300" /></td>
  </tr>
  <tr>
    <td align="center">(a)</td>
    <td align="center">(b)</td>
    <td align="center">(c)</td>
  </tr>
</table>

#### Wave functions

<table>
  <tr>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/.thesis_images/wavefunction/wf2a.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/.thesis_images/wavefunction/wf2b.png" width="300" height="300" /></td>
  </tr>
    <tr>
    <td align="center">(a)</td>
    <td align="center">(b)</td>
  </tr>

</table>

/
### Quantum computing
A size of secular matrix (2n x 2n) is determined by a number of scatteres (n) that creates a system. The computational time for complex systems scales poorly as the number of objects increases.

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
**Email:** migdalekm02@gmail.com
