import numpy as np
import matplotlib.pyplot as plt

from scatterer import scatterer
from struc_const import G_AB
from plots import plot_potential, plot_eigen
from functions import potential, secular, generate_scatterer

#data from article
rA = 0.8
VA = -5
xA = -1
rB = 0.6
VB = -6
RAB = 2
s1 = scatterer(r=rA, V=VA)
s2 = scatterer(rB, VB)
s3 = scatterer(rB, VB)

s = [s1, s2, s3]
tab_x = [xA, xA+RAB, 3]


V , x = potential(scatterers=s, positions=tab_x)
plot_potential(V, x)
plot_eigen(s, tab_x, e1 = -4.5, e2 = -3, n = 100, det = True, values = True)


s, tab_x = generate_scatterer(number=4, min_dist_pot = 0.2, rmin=0.3,  rmax=0.9, Vmin=-2, Vmax=-5, xmin=-4, xmax=4)
V , x = potential(scatterers=s, positions=tab_x)
plot_potential(V, x)
plot_eigen(s, tab_x, e1 = -4.5, e2 = -1, n = 100, det = True, values = True)