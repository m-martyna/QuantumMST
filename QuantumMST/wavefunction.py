from .protein import protein_structure
from .scatterer import scatterer
import numpy as np


# def calc_psi_components(structure):
#     data = {}
#     E = structure.energy
#     data['E'] = E
#     scatterers = structure.scatterers
#     positions = structure.positions
#     alpha = np.sqrt(-E)

#     data['alpha'] = alpha
#     rA, VA = scatterers[0].r, scatterers[0].V
#     rB, VB = scatterers[1].r, scatterers[1].V
#     rC, VC = scatterers[2].r, scatterers[2].V
    
#     RAB = positions[1] - positions[0]
#     RBC = positions[2] - positions[1]
#     RAC = positions[2] - positions[0]
    
#     m0A, m1A = scatterers[0].calc_ml(E, 0), scatterers[0].calc_ml(E, 1)
#     m0B, m1B = scatterers[1].calc_ml(E, 0), scatterers[1].calc_ml(E, 1) 
#     m0C, m1C = scatterers[2].calc_ml(E, 0), scatterers[2].calc_ml(E, 1)
    
#     t0A, t1A = 1/m0A, 1/m1A
#     t0B, t1B = 1/m0B, 1/m1B
#     t0C, t1C = 1/m0C, 1/m1C
    
#     data.update({'t0A': t0A, 't1A': t1A, 't0B': t0B, 't1B': t1B, 't0C': t0C, 't1C': t1C})
    
#     a0A = 1.0 
#     a0B = (np.exp(-alpha * RAB) * (t0A - t1A) + np.exp(-alpha * (RAC + RBC)) * (t0A - t1A) * (t0C - t1C)) / (1.0 - np.exp(-2 * alpha * RBC) * (t0B - t1B) * (t0C - t1C))
#     a0C = (np.exp(-alpha * RAC) * (t0A - t1A) + np.exp(-alpha * (RAB + RBC)) * (t0A - t1A) * (t0B - t1B)) / (1.0 - np.exp(-2 * alpha * RBC) * (t0B - t1B) * (t0C - t1C))
    
#     data.update({'a0A': a0A, 'a0B': a0B, 'a0C': a0C})


#     if E > VA:
#         k_A = np.sqrt(E - VA)
#         data['R0A'], data['R1A'] = np.cos(k_A * rA), np.sin(k_A * rA)
#     else:
#         kappa_A = np.sqrt(VA - E)
#         data['R0A'], data['R1A'] = np.cosh(kappa_A * rA), np.sinh(kappa_A * rA)

#     if E > VB:
#         k_B = np.sqrt(E - VB)
#         data['R0B'], data['R1B'] = np.cos(k_B * rB), np.sin(k_B * rB)
#     else:
#         kappa_B = np.sqrt(VB - E)
#         data['R0B'], data['R1B'] = np.cosh(kappa_B * rB), np.sinh(kappa_B * rB)

#     if E > VC:
#         k_C = np.sqrt(E - VC)
#         data['R0C'], data['R1C'] = np.cos(k_C * rC), np.sin(k_C * rC)
#     else:
#         kappa_C = np.sqrt(VC - E)
#         data['R0C'], data['R1C'] = np.cosh(kappa_C * rC), np.sinh(kappa_C * rC)
    
#     return data


# def calc_psi(x, structure):
#     data = calc_psi_components(structure)
#     E = data['E']
    
#     rA, VA = structure.scatterers[0].r, structure.scatterers[0].V
#     rB, VB = structure.scatterers[1].r, structure.scatterers[1].V
#     rC, VC = structure.scatterers[2].r, structure.scatterers[2].V
#     xA, xB, xC = structure.positions[0], structure.positions[1], structure.positions[2]

#     alpha = data['alpha']
#     rA_x = np.abs(x - xA)
#     rB_x = np.abs(x - xB)
#     rC_x = np.abs(x - xC)
#     sqrt2 = 1/np.sqrt(2)
   
#     if rA_x <= rA:

#         signA = np.sign((xA + rA) - xA)
#         c0A = (data['a0A'] * np.cosh(alpha * rA) + data['a0A'] * data['t0A'] * np.exp(-alpha * rA)) / data['R0A']
#         c1A = signA * (data['a0A'] * np.sinh(alpha * rA) - data['a0A'] * data['t1A'] * np.exp(-alpha * rA)) / data['R1A']
#         c0A *= sqrt2
#         c1A *= sqrt2

#         rA_hat = np.sign(x - xA)
#         R0 = np.cos(np.sqrt(E - VA) * rA_x) if E > VA else np.cosh(np.sqrt(VA - E) * rA_x)
#         R1 = np.sin(np.sqrt(E - VA) * rA_x) if E > VA else np.sinh(np.sqrt(VA - E) * rA_x)
#         return c0A * R0 + c1A * R1 * rA_hat
    

#     elif rB_x <= rB:

#         signB = np.sign((xB - rB) - xB)
#         c0B = (data['a0B'] * np.cosh(alpha * rB) + data['a0B'] * data['t0B'] * np.exp(-alpha * rB)) / data['R0B']
#         c1B = signB * (data['a0B'] * np.sinh(alpha * rB) - data['a0B'] * data['t1B'] * np.exp(-alpha * rB)) / data['R1B']
#         c0B *= sqrt2
#         c1B *= sqrt2

#         rB_hat = np.sign(x - xB)
#         R0 = np.cos(np.sqrt(E - VB) * rB_x) if E > VB else np.cosh(np.sqrt(VB - E) * rB_x)
#         R1 = np.sin(np.sqrt(E - VB) * rB_x) if E > VB else np.sinh(np.sqrt(VB - E) * rB_x)
#         return c0B * R0 + c1B * R1 * rB_hat
    
#     # Region IC (Atom C)
#     elif rC_x <= rC:

#         signC = np.sign((xC - rC) - xC)
#         c0C = (data['a0C'] * np.cosh(alpha * rC) + data['a0C'] * data['t0C'] * np.exp(-alpha * rC)) / data['R0C']
#         c1C = signC * (data['a0C'] * np.sinh(alpha * rC) - data['a0C'] * data['t1C'] * np.exp(-alpha * rC)) / data['R1C']
#         c0C *= sqrt2
#         c1C *= sqrt2

#         rC_hat = np.sign(x - xC)
#         R0 = np.cos(np.sqrt(E - VC) * rC_x) if E > VC else np.cosh(np.sqrt(VC - E) * rC_x)
#         R1 = np.sin(np.sqrt(E - VC) * rC_x) if E > VC else np.sinh(np.sqrt(VC - E) * rC_x)
#         return c0C * R0 + c1C * R1 * rC_hat

#     else:
#         rA_hat = np.sign(x - xA)
#         rB_hat = np.sign(x - xB)
#         rC_hat = np.sign(x - xC)
        
#         psi_A = data['a0A'] * (data['t0A'] - data['t1A'] * rA_hat) * np.exp(-alpha * rA_x)
#         psi_C = data['a0C'] * (data['t0C'] - data['t1C'] * rC_hat) * np.exp(-alpha * rC_x)

#         if x < xB:
#             psi_B = data['a0B'] * (data['t0B'] + data['t1B']) * np.exp(-alpha * rB_x)
#         else:
#             psi_B = data['a0B'] * (data['t0B'] - data['t1B']) * np.exp(-alpha * rB_x)

#         return psi_A + psi_B + psi_C

# def calc_psi_components(structure):
#     data = {}
#     E = structure.energy
#     data['E'] = E
#     scatterers = structure.scatterers
#     positions = structure.positions
#     alpha = np.sqrt(-E)

#     data['alpha'] = alpha
#     rA, VA = scatterers[0].r, scatterers[0].V
#     rB, VB = scatterers[1].r, scatterers[1].V
#     rC, VC = scatterers[2].r, scatterers[2].V
    
#     RAB = positions[1] - positions[0]
#     RBC = positions[2] - positions[1]
#     RAC = positions[2] - positions[0]
#     xA, xB, xC = structure.positions[0], structure.positions[1], structure.positions[2]


#     # Obliczenie m_l oraz t_l
#     m0A, m1A = scatterers[0].calc_ml(E, 0), scatterers[0].calc_ml(E, 1)
#     m0B, m1B = scatterers[1].calc_ml(E, 0), scatterers[1].calc_ml(E, 1) 
#     m0C, m1C = scatterers[2].calc_ml(E, 0), scatterers[2].calc_ml(E, 1)
    
#     t0A, t1A = 1/m0A, 1/m1A
#     t0B, t1B = 1/m0B, 1/m1B
#     t0C, t1C = 1/m0C, 1/m1C
    
#     data.update({'t0A': t0A, 't1A': t1A, 't0B': t0B, 't1B': t1B, 't0C': t0C, 't1C': t1C})
    
#     a0A = 1.0 

#     a0B = (np.exp(-alpha * RAB) * (t0A - t1A) + np.exp(-alpha * (RAC + RBC)) * (t0A - t1A) * (t0C - t1C)) / (1.0 - np.exp(-2 * alpha * RBC) * (t0B - t1B) * (t0C - t1C))
    
#     a0C = (np.exp(-alpha * RAC) * (t0A - t1A) + np.exp(-alpha * (RAB + RBC)) * (t0A - t1A) * (t0B - t1B)) / (1.0 - np.exp(-2 * alpha * RBC) * (t0B - t1B) * (t0C - t1C))
    
#     data.update({'a0A': a0A, 'a0B': a0B, 'a0C': a0C})

#     if E > VA:
#         k_A = np.sqrt(E - VA)
#         R0A, R1A = np.cos(k_A * rA), np.sin(k_A * rA)
#     else:
#         kappa_A = np.sqrt(VA - E)
#         R0A, R1A = np.cosh(kappa_A * rA), np.sinh(kappa_A * rA)

#     # Atom B
#     if E > VB:
#         k_B = np.sqrt(E - VB)
#         R0B, R1B = np.cos(k_B * rB), np.sin(k_B * rB)
#     else:
#         kappa_B = np.sqrt(VB - E)
#         R0B, R1B = np.cosh(kappa_B * rB), np.sinh(kappa_B * rB)

#     # Atom C
#     if E > VC:
#         k_C = np.sqrt(E - VC)
#         R0C, R1C = np.cos(k_C * rC), np.sin(k_C * rC)
#     else:
#         kappa_C = np.sqrt(VC - E)
#         R0C, R1C = np.cosh(kappa_C * rC), np.sinh(kappa_C * rC)


#     c0A = a0A * np.cosh(alpha * rA) + a0A * t0A * np.exp(-alpha * rA)

#     c1A = a0A * np.sinh(alpha * rA) - a0A * t1A * np.exp(-alpha * rA)

#     data['c0A'], data['c1A'] = c0A / R0A, c1A / R1A



#     # Atom B (Środkowy dołek - odwrócenie fazy fali antysymetrycznej, równanie 28)

#     c0B = a0B * np.cosh(alpha * rB) + a0B * t0B * np.exp(-alpha * rB)

#     c1B = -a0B * np.sinh(alpha * rB) + a0B * t1B * np.exp(-alpha * rB) # POPRAWIONE ZNAKI

#     data['c0B'], data['c1B'] = c0B / R0B, c1B / R1B



#     # Atom C (Kolejny dołek po prawej stronie - analogicznie do B zachowuje odwrócenie fazy)

#     c0C = a0C * np.cosh(alpha * rC) + a0C * t0C * np.exp(-alpha * rC)

#     c1C = -a0C * np.sinh(alpha * rC) - a0C * t1C * np.exp(-alpha * rC) # POPRAWIONE ZNAKI

#     data['c0C'], data['c1C'] = c0C / R0C, c1C / R1C 
    
#     return data


# def calc_psi(x, structure):
#     data = calc_psi_components(structure)
#     E = data['E']
    
#     rA, VA = structure.scatterers[0].r, structure.scatterers[0].V
#     rB, VB = structure.scatterers[1].r, structure.scatterers[1].V
#     rC, VC = structure.scatterers[2].r, structure.scatterers[2].V
#     xA, xB, xC = structure.positions[0], structure.positions[1], structure.positions[2]

#     alpha = data['alpha']
#     rA_x = np.abs(x - xA)
#     rB_x = np.abs(x - xB)
#     rC_x = np.abs(x - xC)

#     # Region IA
#     if rA_x <= rA:
#         rA_hat = np.sign(x - xA)
#         R0 = np.cos(np.sqrt(E - VA) * rA_x) if E > VA else np.cosh(np.sqrt(VA - E) * rA_x)
#         R1 = np.sin(np.sqrt(E - VA) * rA_x) if E > VA else np.sinh(np.sqrt(VA - E) * rA_x)
#         return data['c0A'] * R0 + data['c1A'] * R1 * rA_hat
    
#     # Region IB
#     elif rB_x <= rB:
#         rB_hat = np.sign(x - xB)
#         R0 = np.cos(np.sqrt(E - VB) * rB_x) if E > VB else np.cosh(np.sqrt(VB - E) * rB_x)
#         R1 = np.sin(np.sqrt(E - VB) * rB_x) if E > VB else np.sinh(np.sqrt(VB - E) * rB_x)
#         return data['c0B'] * R0 + data['c1B'] * R1 * rB_hat
    
#     # Region IC
#     elif rC_x <= rC:
#         rC_hat = np.sign(x - xC)
#         R0 = np.cos(np.sqrt(E - VC) * rC_x) if E > VC else np.cosh(np.sqrt(VC - E) * rC_x)
#         R1 = np.sin(np.sqrt(E - VC) * rC_x) if E > VC else np.sinh(np.sqrt(VC - E) * rC_x)
#         return data['c0C'] * R0 + data['c1C'] * R1 * rC_hat


#     else:
#         rA_hat = np.sign(x - xA)
#         rB_hat = np.sign(x - xB)
#         rC_hat = np.sign(x - xC)
        
#         psi_A = data['a0A'] * (data['t0A'] - data['t1A'] * rA_hat) * np.exp(-alpha * rA_x)
        
#         psi_C = data['a0C'] * (data['t0C'] - data['t1C'] * rC_hat) * np.exp(-alpha * rC_x)

#         if x < xB:
#             psi_B = data['a0B'] * (data['t0B'] + data['t1B']) * np.exp(-alpha * rB_x)
#         else:
#             psi_B = data['a0B'] * (data['t0B'] - data['t1B']) * np.exp(-alpha * rB_x)

#         return psi_A + psi_B + psi_C













# def calc_psi_components(structure: protein_structure):
    
#     data = {}
#     E = structure.energy
#     data['E'] = E
#     scatterers = structure.scatterers
#     positions = structure.positions
#     alpha = np.sqrt(-E)

#     data['alpha'] = alpha
#     rA, VA = scatterers[0].r, scatterers[0].V
#     rB, VB = scatterers[1].r, scatterers[1].V
#     rC, VC = scatterers[2].r, scatterers[2].V
#     RAB = positions[1] - positions[0]
#     RBC = positions[2] - positions[1]
#     RAC = positions[2] - positions[0]
#     RBC = positions[2] - positions[1]


#     m0A, m1A = scatterers[0].calc_ml(E, 0), scatterers[0].calc_ml(E, 1)
#     m0B, m1B = scatterers[1].calc_ml(E, 0), scatterers[1].calc_ml(E, 1) 
#     m0C, m1C = scatterers[2].calc_ml(E, 0), scatterers[2].calc_ml(E, 1)
#     t0A = 1/m0A
#     t1A = 1/m1A
#     t0B = 1/m0B
#     t1B = 1/m1B
#     t0C = 1/m0C
#     t1C = 1/m1C
    
#     data['t0A'] = t0A
#     data['t1A'] = t1A
#     data['t0B'] = t0B
#     data['t1B'] = t1B
#     data['t0C'] = t0C
#     data['t1C'] = t1C
    
#     a0A = 1 #assumption
#     #a0B = np.exp(-alpha*RAB) * (t0A-t1A)*a0A #from (24b) and (11)

#     AB = np.exp(-alpha * RAB) * (t0A - t1A)
#     AC = np.exp(-alpha * RAC) * (t0A - t1A)  
#     BC = np.exp(-alpha * RBC) * (t0B - t1B)
#     CB = np.exp(-alpha * RBC) * (t0C + t1C)

#     a0B = (AB + CB * AC) / (1.0 - (BC * CB))
#     a0C = (AC + BC * AB) / (1.0 - (BC * CB))
    
    
#     data['a0A'] = a0A
#     data['a0B'] = a0B
#     data['a0C'] = a0C

#     #Region IA
#     if(E>VA):
#         k_A = np.sqrt(E-VA)
#         R0A = np.cos(k_A*rA)
#         R1A = np.sin(k_A*rA)

#         data['k_A'] = k_A
#     else:
#         kappa_A = np.sqrt(VA-E)
#         R0A = np.cosh(kappa_A*rA)
#         R1A = np.sinh(kappa_A*rA)

#         data['kappa_A'] = kappa_A
        
#     #Region IB
#     if(E>VB):
#         k_B = np.sqrt(E-VB)
#         R0B = np.cos(k_B*rB)
#         R1B = np.sin(k_B*rB)

#         data['k_B'] = k_B
#     else:
#         kappa_B = np.sqrt(VB-E)
#         R0B = np.cosh(kappa_B*rB)
#         R1B = np.sinh(kappa_B*rB)

#         data['kappa_B'] = kappa_B

#         #Region IC
#     if(E>VC):
#         k_C = np.sqrt(E-VC)
#         R0C = np.cos(k_C*rC)
#         R1C = np.sin(k_C*rC)

#         data['k_C'] = k_C
#     else:
#         kappa_C = np.sqrt(VC-E)
#         R0C = np.cosh(kappa_C*rC)
#         R1C = np.sinh(kappa_C*rC)

#         data['kappa_C'] = kappa_C

#     c0A = a0A*np.cosh(alpha*rA) + a0A*t0A*np.exp(-alpha*rA) #from (3) and (27)
#     c1A = a0A*np.sinh(alpha*rA) - a0A*t1A*np.exp(-alpha*rA)
    
#     data['c0A'] = c0A/R0A
#     data['c1A'] = c1A/R1A

#     c0B = a0B*np.cosh(alpha*rB) + a0B*t0B*np.exp(-alpha*rB)
#     c1B = -a0B*np.sinh(alpha*rB) + a0B*t1B* np.exp(-alpha*rB)

#     data['c0B'] = c0B/R0B
#     data['c1B'] = c1B/R1B

#     c0C = a0C*np.cosh(alpha*rC) + a0C*t0C*np.exp(-alpha*rC)
#     c1C = a0C*np.sinh(alpha*rC) + a0C*t1C* np.exp(-alpha*rC)

#     data['c0C'] = c0C/R0C
#     data['c1C'] = c1C/R1C
    
#     return data


# def calc_psi(x, structure: protein_structure):

#     data = calc_psi_components(structure)
    
#     E = data['E']
#     rA, VA = structure.scatterers[0].r, structure.scatterers[0].V
#     rB, VB = structure.scatterers[1].r, structure.scatterers[1].V
#     rC, VC = structure.scatterers[2].r, structure.scatterers[2].V
#     xA, xB, xC = structure.positions[0], structure.positions[1], structure.positions[2]

#     alpha = data['alpha']
#     rA_x = np.abs(x-xA)
#     rB_x = np.abs(x-xB)
#     rC_x = np.abs(x-xC)

#     #Region IA
#     if rA_x <= rA:
#         rA_hat = np.sign(x-xA)
#         if E > VA:
#             R0 = np.cos(data['k_A']*rA_x)
#             R1 = np.sin(data['k_A']*rA_x)
#         else:
#             R0 = np.cosh(data['kappa_A']*rA_x)
#             R1 = np.sinh(data['kappa_A']*rA_x)
#         return data['c0A']*R0 + data['c1A']*R1* rA_hat
    
#     #Region IB
#     elif rB_x <= rB:
#         rB_hat = np.sign(x-xB)
#         if E > VB:
#             R0 = np.cos(data['k_B']*rB_x)
#             R1 = np.sin(data['k_B']*rB_x)
#         else:
#             R0 = np.cosh(data['kappa_B']*rB_x)
#             R1 = np.sinh(data['kappa_B']*rB_x)
#         return data['c0B']*R0 + data['c1B']*R1* rB_hat
    
#     #Region IC
#     elif rC_x <= rC:
#         rC_hat = np.sign(x-xC)
#         if E > VC:
#             R0 = np.cos(data['k_C']*rC_x)
#             R1 = np.sin(data['k_C']*rC_x)
#         else:
#             R0 = np.cosh(data['kappa_C']*rC_x)
#             R1 = np.sinh(data['kappa_C']*rC_x)
#         return data['c0C']*R0 + data['c1C']*R1* rC_hat

#     #Region II
#     else:
#         rA_hat = np.sign(x - xA)
#         rB_hat = np.sign(x - xB)
#         rC_hat = np.sign(x - xC)
        
#         psi_A = data['a0A'] * (data['t0A'] - data['t1A'] * rA_hat) * np.exp(-alpha * rA_x)
#         psi_B = data['a0B'] * (data['t0B'] - data['t1B'] * rB_hat) * np.exp(-alpha * rB_x)
#         psi_C = data['a0C'] * (data['t0C'] - data['t1C'] * rC_hat) * np.exp(-alpha * rC_x)
        
#         return psi_A + psi_B + psi_C









def calc_psi_components(structure: protein_structure):
    
    data = {}
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


