import numpy as np
import matplotlib.pyplot as plt

from scatterer import scatterer
from struc_const import G_AB
from plots import plot_potential, plot_eigen
from functions import potential, secular, generate_scatteters, generate_positions

#data from article
rA = 0.8
VA = -5
xA = -1
rB = 0.6
VB = -6
RAB = 2
s1 = scatterer(r = rA, V = VA)
s2 = scatterer(rB, VB)
s = [s1, s2]
tab_x = [xA, xA+RAB]

# s = generate_scatteters(4, -4, -2, 0.2, 2)
# tab_x = generate_positions(s, -40, 40)

V , x = potential(scatterers = s, positions = tab_x)
plot_potential(V, x)

plot_eigen(s, tab_x, e1 = -4.5, e2 = -3, n = 100, determinans = True, values = True)
