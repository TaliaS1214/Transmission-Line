# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/Users/SeanTalia/.spyder2/.temp.py
"""

import numpy as np

def main():
    
    m = 4 # Numer of Segments
    N = 2*m + 1 # Dimension of System     
    
    R = 1.0    
    L = 1.0
    G = 0.5
    C = 0.3
    
    omega = 2.0     # Voltage Frequency
    V = 5.0         # Voltage Amplitude
    
    Z = R + omega*L*(1.0j)  # Line series impedance
    Y = G + omega*C*(1.0j)  # Line shunt admittance
    
    Z_L = 5.0 + 3.0j  # Load impedance

    # Below vectors are of the form [V_1, V_2, V_3, V_4, I_1, I_2, I_3, I_4, I_L]
    
    E_1 = [ 1.0,     0,     0,     0,     Z,     0,     0,      0,     0]
    E_2 = [-1.0,   1.0,     0,     0,     0,     Z,     0,      0,     0]
    E_3 = [   0,  -1.0,   1.0,     0,     0,     0,     Z,      0,     0]
    E_4 = [   0,     0,  -1.0,   1.0,     0,     0,     0,      Z,     0]
    
    E_5 = [   Y,     0,     0,     0,  -1.0,   1.0,     0,      0,     0]    
    E_6 = [   0,     Y,     0,     0,     0,  -1.0,   1.0,      0,     0]
    E_7 = [   0,     0,     Y,     0,     0,     0,  -1.0,    1.0,     0]
    E_8 = [   0,     0,     0,     Y,     0,     0,     0,   -1.0,   1.0]    
    
    E_9 = [   0,     0,     0,     0,     Z,     Z,     Z,     Z,    Z_L]
    
    A = np.matrix([E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8, E_9])
    
    b = np.zeros(N)    
    b[0] = V
    b[N-1] = V
    
    sol = np.linalg.solve(A, b)
    
    for i in range(N):
        if i == 0:
            print "\n Voltage Constants: "
        if i == np.floor(N/2):
            print "\n Current through Branches: "
        if i == N-1:
            print "\n Current through load: "
        print "Constant: ", sol[i], "Magnitude: ", abs(sol[i])

main()