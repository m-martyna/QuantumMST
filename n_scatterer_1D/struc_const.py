import numpy as np

class G_AB:
    def __init__(self, energy, positions):
        self.E = energy
        self.pos = positions

    def calc_g(self, xA, xB):
        XAB = xB - xA
        RAB = np.abs(XAB)

        g = np.zeros((2, 2))
        phase = np.exp(1j * np.sqrt(complex(self.E)) * RAB)

        g = np.array([[phase, 1j * np.sign(XAB) * phase],
                      [-1j * np.sign(XAB) * phase, phase]],
                      dtype=complex)

        return g
    
    def matrix_g(self):
        range_x = len(self.pos)
        matrix_g = np.zeros((range_x , range_x, 2, 2 ), dtype = complex)
        for i in range (0, range_x):
            for j in range (0, range_x):
                x1 = self.pos[i]
                x2 = self.pos[j]
                matrix_g[i, j] = self.calc_g(x1, x2)

        return matrix_g

