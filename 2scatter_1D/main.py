import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from scatters import two_scatters, secular, calc_psi_components, calc_psi

S2 = two_scatters(0.8, -5, -1, 0.6, -6, 2, 4, analytical = True)
V, array_x = S2.potential()

def calc_energies(S2):
    eigenenergies = []

    E_range = np.arange(min(S2.VA, S2.VB), 0., 0.3)
    VA_VB = (np.abs(E_range - S2.VA) >= 0.2) & (np.abs(E_range - S2.VB) >= 0.2)
    E_range = E_range[VA_VB]

    e1 = fsolve(secular, E_range[0], args=(S2.rA, S2.VA, S2.rB, S2.VB, S2.RAB))[0]
    eigenenergies.append(e1)

    for e in E_range[1:]:
        e2 = fsolve(secular, e, args=(S2.rA, S2.VA, S2.rB, S2.VB, S2.RAB))[0]
        if not (np.isclose(e1,e2, atol = 0.1)): 
            eigenenergies.append(e2)
            e1 = e2
    return eigenenergies

energies = calc_energies(S2)

print("Calculated energies:")
print(f"E = {energies[0]:.9f}")
print(f"E = {energies[1]:.9f}")
print("\nEnergies FIG 2.:")
print("E = -3.949901757")
print("E = -3.22848704")


components_e1 = calc_psi_components(S2,energies[0])
components_e2 = calc_psi_components(S2, energies[1])

psi_e1 = np.array([calc_psi(S2, x, components_e1) for x in array_x])
psi_e2 = np.array([calc_psi(S2, x, components_e2) for x in array_x])

import scipy.integrate 
integral_e1 = scipy.integrate.simpson((np.abs(psi_e1)**2), x = array_x)
integral_e2 = scipy.integrate.simpson((np.abs(psi_e2)**2), x = array_x)
norm_psi_e1 = psi_e1/np.sqrt(integral_e1)
norm_psi_e2 = psi_e2/np.sqrt(integral_e2)
a = np.max(norm_psi_e1)*1.1

plt.figure()
plt.plot(array_x, norm_psi_e1, color='black')
plt.axhline(0, ls='--', color='black')
plt.fill_between(array_x, -a, a, where=V<0, color='red', alpha=0.1)
plt.ylim(-a, a)
plt.title(f"E = {energies[0]}")
plt.xlabel("x")
plt.ylabel(f"$\Psi$")
plt.savefig("2scatter_1D/plots/Psi1.png")

plt.figure()
plt.plot(array_x, norm_psi_e2, color='black') 
plt.axhline(0, ls='--',  color='black')
plt.fill_between(array_x, -a, a, where=V<0, color='red', alpha=0.1)
plt.ylim(-a, a)
plt.title(f"E = {energies[1]}")
plt.xlabel("x")
plt.ylabel(f"$\Psi$")
plt.savefig("2scatter_1D/plots/Psi2.png")