import numpy as np
import matplotlib.pyplot as plt

def calc_R(E, V, r, l, analitical = True):

    if(analitical):
        if(E>V):
            k = np.sqrt(E-V)
            if(l==0): R = np.cos(k*r)
            if(l==1): R = np.sin(k*r)
        else:
            kappa = np.sqrt(V-E)
            if(l==0): R = np.cosh(kappa*r)
            if(l==1): R = np.sinh(kappa*r)

    return R

def calc_Rderivative(E, V, r, l, analitical = True):

    if(analitical):
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
    

class two_scatters:
    def __init__(self, rA, VA, xA, rB, VB, RAB, x_range, analytical = True):
        self.rA = rA
        self.VA = VA
        self.xA = xA
        self.rB = rB
        self.VB = VB
        self.RAB = RAB
        self.x_range = x_range
        self.x_range
        self.xB = self.xA+self.RAB

    def potential(self, plot = True):
        array_x = np.linspace(-self.x_range, self.x_range, 100)
        V = []

        for x in array_x:
            if (x>=self.xA-self.rA) and (x<=self.xA+self.rA):
                V.append(self.VA)
            elif (x>=self.xB-self.rB) and (x<=self.xB+self.rB):
                V.append(self.VB)
            else:
                V.append(0)

        plt.plot(array_x, V, color = "black")
        plt.title("V(x)")
        plt.xlabel("x")
        plt.ylabel("V")
        plt.xlim(-self.x_range, self.x_range)
        plt.ylim(self.VB*1.5, np.abs(self.VB)*0.75)
        if(plot):
            plt.savefig("2scatter_1D/plots/scatters_illustration.png")

        return np.array(V), array_x

def gamma(E, V, r):

    R0 = calc_R(E, V, r, 0, analitical = True)
    R0_p = calc_Rderivative(E, V, r, 0, analitical = True)
    
    R1 = calc_R(E, V, r, 1, analitical = True)
    R1_p = calc_Rderivative(E, V, r, 1, analitical = True)
    
    
    return R0_p/R0, R1_p/R1

def m_matrix(E, gamma_0, gamma_1, r0):

    if E >= 0:
        return np.nan, np.nan
    
    k = np.sqrt(E + 0j) 
    z = k * r0 
    
    den0 = z * j_prime_1d(0, z) - r0 * gamma_0 * j_1d(0, z)
    num0 = z * h_prime_1d(0, z) - r0 * gamma_0 * h_1d(0, z) 

    den1 = z * j_prime_1d(1, z) - r0 * gamma_1 * j_1d(1, z)
    num1 = z * h_prime_1d(1, z) - r0 * gamma_1 * h_1d(1, z)

    return -num0/den0, -num1/den1


def secular(E, r_A, V_A, r_B, V_B, R_AB):

    alpha = np.sqrt(-E)

    gamma_0A, gamma_1A = gamma(E, V_A, r_A)
    gamma_0B, gamma_1B = gamma(E, V_B, r_B)

    m0A, m1A = m_matrix(E, gamma_0A, gamma_1A, r_A)
    m0B, m1B = m_matrix(E, gamma_0B, gamma_1B, r_B)

    term1 = np.exp(-2*alpha*R_AB) * (m1A-m0A) * (m1B-m0B)
    term2 = m0A*m1A*m0B*m1B

    result = term1-term2

    return np.real(result)

def calc_psi_components(S, E):
    rA = S.rA
    VA = S.VA
    rB = S.rB
    VB = S.VB
    RAB = S.RAB
    
    data = {}

    data['E'] = E

    alpha = np.sqrt(-E)

    data['alpha'] = alpha

    g0A, g1A = gamma(E, VA, rA)
    g0B, g1B = gamma(E, VB, rB)    
    m0A, m1A = m_matrix(E, g0A, g1A, rA)
    m0B, m1B = m_matrix(E, g0B, g1B, rB)   
    t0A = 1/m0A
    t1A = 1/m1A
    t0B = 1/m0B
    t1B = 1/m1B
    
    data['t0A'] = t0A
    data['t1A'] = t1A
    data['t0B'] = t0B
    data['t1B'] = t1B
    
    a0A = 1 #assumption
    a0B = np.exp(-alpha*RAB) * (t0A-t1A)*a0A #from (24b) and (11)
    
    data['a0A'] = a0A
    data['a0B'] = a0B

    #Region IA
    if(E>VA):
        k_A = np.sqrt(E-VA)
        R0A = np.cos(k_A*rA)
        R1A = np.sin(k_A*rA)

        data['k_A'] = k_A
    else:
        kappa_A = np.sqrt(VA-E)
        R0A = np.cosh(kappa_A*rA)
        R1A = np.sinh(kappa_A*rA)

        data['kappa_A'] = kappa_A
        
    #Region IB
    if(E>VB):
        k_B = np.sqrt(E-VB)
        R0B = np.cos(k_B*rB)
        R1B = np.sin(k_B*rB)

        data['k_B'] = k_B
    else:
        kappa_B = np.sqrt(VB-E)
        R0B = np.cosh(kappa_B*rB)
        R1B = np.sinh(kappa_B*rB)

        data['kappa_B'] = kappa_B

    c0A = a0A*np.cosh(alpha*rA) + a0A*t0A*np.exp(-alpha*rA) #from (3) and (27)
    c1A = a0A*np.sinh(alpha*rA) - a0A*t1A*np.exp(-alpha*rA)
    
    data['c0A'] = c0A/R0A
    data['c1A'] = c1A/R1A

    c0B = a0B*np.cosh(alpha*rB) + a0B*t0B*np.exp(-alpha*rB)
    c1B = -a0B*np.sinh(alpha*rB) + a0B*t1B* np.exp(-alpha*rB)

    data['c0B'] = c0B/R0B
    data['c1B'] = c1B/R1B
    
    return data

def calc_psi(S, x, data):
    rA = S.rA
    VA = S.VA
    xA = S.xA
    xB = S.xB
    rB = S.rB
    VB = S.VB
    RAB = S.RAB

    E = data['E']
    alpha = data['alpha']
    rA_x = np.abs(x-xA)
    rB_x = np.abs(x-xB)

    #Region IA
    if rA_x <= rA:
        rA_hat = np.sign(x-xA)
        if E > VA:
            R0 = np.cos(data['k_A']*rA_x)
            R1 = np.sin(data['k_A']*rA_x)
        else:
            R0 = np.cosh(data['kappa_A']*rA_x)
            R1 = np.sinh(data['kappa_A']*rA_x)
        return data['c0A']*R0 + data['c1A']*R1* rA_hat
    
    #Region IB
    elif rB_x <= rB:
        rB_hat = np.sign(x-xB)
        if E > VB:
            R0 = np.cos(data['k_B']*rB_x)
            R1 = np.sin(data['k_B']*rB_x)
        else:
            R0 = np.cosh(data['kappa_B']*rB_x)
            R1 = np.sinh(data['kappa_B']*rB_x)
        return data['c0B']*R0 + data['c1B']*R1* rB_hat

    #Region II
    else:
        rA_hat = np.sign(x-xA)
        rB_hat = np.sign(x-xB)
        
        psi_A = data['a0A']*(data['t0A'] - data['t1A']*rA_hat) * np.exp(-alpha*rA_x)
        psi_B = data['a0B']*(data['t0B'] + data['t1B']*rB_hat) * np.exp(-alpha*rB_x)
        
        return psi_A + psi_B

 