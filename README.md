# **Application of quantum computing to electronic structure calculation using Green’s function method.**

This project aims to explore multiple scattering theory and its potential applications within the field of quantum computing.


## Table of contents
* [Butler’s approach to MST](#butler)
* [Validation of results from article](#validation)



## Butler’s Approach to MST
The core of calculations is based on W. H. Butler’s article _Validity and accuracy of multiple-scattering theory_ (https://doi.org/10.1103/PhysRevB.42.1518). The primary objective is to determine the solution to the wave equation for a muffin-tin potential. This potenitial in one dimension is a sum of symmetric, disjoint potentials

$$V(x) = \sum_{n} v_n(x - X_n).$$

There are two distinct regions (I and II) with non-vanishing and zero potential. The wave functions in these regions are given by

$$\Psi_{In} = \sum_{l=0,1} c_l^n R_l^n(r_n) Y_l(\hat{r}_n),$$

$$\Psi_{II} = \sum_{n=A,B} \sum_{l=0,1} b_l^n h_l(E^{1/2} r_n) Y_l(\hat{r}_n).$$

<table>
  <tr>
  <td><img src="https://github.com/m-martyna/Praca_magisterska/blob/main/examples/plots/scatters_illustration.png" width="400" height="400" /></td>
  </tr>
  <tr>
    <td align="center"></td>
  </tr>
</table>

## Validation of results from article
As an initial benchmark I applied the theory to a one-dimensional two-scatterer problem following the approach described in Butler’s article (https://doi.org/10.1103/PhysRevB.42.1518).

<table>
  <tr>
  <td><img src="https://github.com/m-martyna/Praca_magisterska/blob/main/examples/plots/Psi1.png" width="400" height="400" /></td>
  <td><img src="https://github.com/m-martyna/Praca_magisterska/blob/main/examples/plots/Psi2.png" width="400" height="400" /></td>
  </tr>
  <tr>
    <td align="center">(a)</td>
    <td align="center">(b)</td>
  </tr>
</table>

/home/ugi/Praca_magisterska/examples/plots/Psi1.png