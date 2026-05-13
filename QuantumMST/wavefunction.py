from .protein import protein_structure
from .scatterer import scatterer
import numpy as np


def calc_psi_components(structure: protein_structure):
    data = {}

    print(type(structure))
    E = structure.energy
    data['E'] = E
    scatterers = structure.scatterers
    positions = structure.positions
    alpha = np.sqrt(-E)

    data['alpha'] = alpha
    rA, VA = scatterers[0].r, scatterers[0].V
    rB, VB = scatterers[1].r, scatterers[1].V
    RAB = positions[1] - positions[0]


    m0A, m1A = scatterers[0].calc_ml(E, 0), scatterers[0].calc_ml(E, 1)
    m0B, m1B = scatterers[1].calc_ml(E, 0), scatterers[1].calc_ml(E, 1) 
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


def calc_psi(x, structure: protein_structure):

    data = calc_psi_components(structure)
    
    E = data['E']
    rA, VA = structure.scatterers[0].r, structure.scatterers[0].V
    rB, VB = structure.scatterers[1].r, structure.scatterers[1].V
    xA, xB = structure.positions[0], structure.positions[1]

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