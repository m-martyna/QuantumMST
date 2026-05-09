# **Quantum computing in MST**

This project aims to explore multiple scattering theory and its potential applications within the field of quantum computing. This research was conducted as part of my Master’s degree titled **_Application of quantum computing to electronic structure calculation using Green’s function method_**, which provides more complex analysis of the topics covered in this repository.

## Requirements
The code was developed and tested using the following versions of:
* **Python:** `3.12.2`
* **Qiskit:** `2.2.3`

### Installation
Ensure you have Python installed ([python.org](https://www.python.org/)). The required quantum computing framework can be set up via `pip`:
```bash
pip install qiskit==2.2.3
```
## Overview
The core of calculations is based on W. H. Butler’s article
> **_Validity and accuracy of multiple-scattering theory_**  
> https://doi.org/10.1103/PhysRevB.41.2684.   

The scattering problem for a 1D system with disjoint, symmetric potentials is equivalent to solving the secular matrix eigenvalue problem. If a given energy results in a zero determinant for the system, it is considered an allowed energy level. This matrix is constructed using a list of **`scatterer`** objects initialized with the dimensions of each well, an array of the center positions of scatterers and a chosen energy.  

To verify the implementation of the matrix, the analytical results from the article were compared with numerical results obtained using the NumPy library. The potential (a), the determinant  (b) as well as the eigenvalues (c) for various energy values are plotted below.

<table>
  <tr>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/Butler/potential.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/Butler/determinant.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/Butler/eigen_values.png" width="300" height="300" /></td>
  </tr>
  <tr>
    <td align="center">(a)</td>
    <td align="center">(b)</td>
    <td align="center">(c)</td>
  </tr>
</table>

The wave functions plotted for the determined energies (a)(b) perfectly match those provided in the article (c)(d).

<table>
  <tr>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/Butler/Psi1.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/Butler/Psi2.png" width="300" height="300" /></td>
  </tr>
  <tr>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/Butler/Psi1_Butler.png" width="300" height="300" /></td>
  <td><img src="https://github.com/m-martyna/QuantumMST/blob/main/results/images/Butler/Psi2_Butler.png" width="300" height="300" /></td>
  </tr>
  <tr>
    <td align="center">(a)</td>
    <td align="center">(b)</td>
  </tr>
</table>




## Contact & attribution
**Author:** Martyna Migdałek    
**Supervised by**:  
* **Tomasz Stopa, PhD**
* **Jakub Haberko, PhD**  

**Institution:** AGH University of Krakow, Faculty of Physics and Applied Computer Science  
**Email:** migdalekm@student.agh.edu.pl
