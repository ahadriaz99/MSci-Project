# -*- coding: utf-8 -*-
""""
Testing script for identifying incompressible liquid phases at half-fillings


25/10/21

by Marcell Dorian Kovacs and Ahad Riaz
"""

from Hamiltonian import Hamiltonian
import numpy as np
import math
import matplotlib.pyplot as plt
import copy

from Fock_vector import fock_vector
import Ryser_Algorithm as ryser
import Sphere_Hamiltonian as sphere
import Disc_Hamiltonian as disc 
import config as configs
from numpy import linalg as la

N0 = 10
S = 4
M = 5
mu = N0/(2*S)
#N_range = np.array([N0-1,N0, N0+1])
e_grounds = []      
e_values = []
eprime_grounds = [] # For sign problems
L_range = np.linspace(0,M,M+1)
gap = []

eplus = 0
eminus = 0 
e = 0

'''for N in range(N0):
    if N == 0 or N == 1:
        continue
    else:'''
for N in range(N0+1):
    print(N)
    N_range = np.array([N-1,N, N+1])
    if N == 0 or N == 1:
        gap.append(0)
    else:
        for i in N_range:
            H = disc.disc_Hamiltonian(N = i, M = M, L = 6)
            H.generate_basis()
            H.construct_Hamiltonian()
            evalues,evecs = H.diagonalise()
            if i == N + 1:
                eplus = H.e_ground/i
            if i == N - 1:
                eminus = H.e_ground/i
            if i == N:
                e = H.e_ground/i
                                        
        gap.append(N*(eplus + eminus - 2*e))

plt.figure()
plt.grid()
plt.plot(range(N0+1),gap, 'x')
plt.show()
        