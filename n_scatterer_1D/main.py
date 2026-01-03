import numpy as np
import matplotlib.pyplot as plt
import time

from scatterer import scatterer
from struc_const import G_AB
from functions import potential, secular, eigen_val, generate_scatterer, load_scatterer, interpolation, generate_even_scatterer
#! data from Butler's article
rA = 0.8
VA = -5
xA = -1
rB = 0.6
VB = -6
RAB = 2


#! initializate scatterers and their coordinates
s1 = scatterer(r=rA, V=VA)
s2 = scatterer(rB, VB)

s3 = scatterer(rB, VB)
s = [s1, s2]
tab_x = [xA, xA+RAB]

# s, tab_x = generate_scatterer(number=5, min_dist_pot = 0.5,
#                                rmin=0.5,  rmax=1.9,
#                                 Vmin=-2, Vmax=-5,
#                                 xmin=-10, xmax=10, save = True)

# s, tab_x = load_scatterer('n_scatterer_1D/results/ex.txt')
# s, tab_x = load_scatterer('n_scatterer_1D/results/random_scatterer_5.txt')


#! calculate and plot potential, eigenvalues, and determinant
nE = 50
V = potential(scatterers=s, positions=tab_x)
eigen_val(s, tab_x, e1 = -5.5, e2 = -2.5, n = nE,
           plot_det = False, plot_values = True)

#! interpolate eigenvalues to reduce computational cost
interpolation(nE)


#! estimate computational complexity
# array_time = []
# array_n = range(2, 100) #* range of the number of scatterers
# E = 3.5
# for n in array_n:
#     print(n)

#     # s, tab_x = generate_scatterer(number=n, min_dist_pot = 0.5,
#     #                            rmin=1.25,  rmax=1.5,
#     #                             Vmin=-5, Vmax=-5,
#     #                             xmin=-2.25*n, xmax=2.25*n, save = False)
#     s, tab_x = generate_even_scatterer(n)
    
#     # V = potential(scatterers=s, positions=tab_x)
#     str_constans = G_AB(energy = E , positions = tab_x)
    
#     start_time = time.perf_counter()
#     eigen_values = secular(s, str_constans, only_eigen_val=True)
#     end_time = time.perf_counter()
#     array_time.append(end_time-start_time)

# plt.figure()
# plt.plot(array_n, array_time, 'o-', color = 'black', markersize=2)
# plt.xlabel("number of scatterers")
# plt.ylabel("time[s]")

# plt.show()