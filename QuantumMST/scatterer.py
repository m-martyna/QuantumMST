import numpy as np

def calc_R(E, V, r, l):
    if(E>V):
        k = np.sqrt(E-V)
        if(l==0): R = np.cos(k*r)
        if(l==1): R = np.sin(k*r)
    else:
        kappa = np.sqrt(V-E)
        if(l==0): R = np.cosh(kappa*r)
        if(l==1): R = np.sinh(kappa*r)

    return R

def calc_Rderivative(E, V, r, l):
    if(E>V):
        k = np.sqrt(E-V)
        if(l==0): R_prim = -k*np.sin(k*r)
        if(l==1): R_prim = k*np.cos(k*r)
    else:
        kappa = np.sqrt(V-E)
        if(l==0): R_prim = kappa*np.sinh(kappa*r)
        if(l==1): R_prim = kappa*np.cosh(kappa*r)

    return R_prim

def j_1d(l, z):
    return np.cos(z) if (l==0) else np.sin(z)

def j_prime_1d(l, z):
    return -np.sin(z) if (l==0) else np.cos(z)

def h_1d(l, z):
    return (np.cos(z) + 1j * np.sin(z)) if l == 0 else (np.sin(z) - 1j * np.cos(z))

def h_prime_1d(l, z):
    return (-np.sin(z) + 1j * np.cos(z)) if l == 0 else (np.cos(z) + 1j * np.sin(z))

class scatterer:
    def __init__(self, r, V):
        self.r = r
        self.V = V

    def calc_gamma(self, E, l):
        R = calc_R(E, self.V, self.r, l)
        R_p = calc_Rderivative(E, self.V, self.r, l)
        
        
        epsilon = 1e-14
        if abs(R) < epsilon:
            if abs(R_p) < epsilon:
                if self.r != 0:
                    return 1.0 / self.r
                else:
                    return 0.0
                
        return R_p/R
    
    def calc_ml(self, E, l):        
        k = np.sqrt(E + 0j) 
        z = k * self.r 
        gamma = self.calc_gamma(E, l)
        denominator = z * j_prime_1d(l, z) - self.r * gamma * j_1d(l, z)
        numerator = z * h_prime_1d(l, z) - self.r * gamma * h_1d(l, z) 
        return -numerator/denominator
    
    def calc_m(self, E):
        m = np.array([[self.calc_ml(E, 0), 0],
                      [0, self.calc_ml(E, 1)]],
                      dtype=complex)
        
        return m