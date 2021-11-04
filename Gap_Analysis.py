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
import config as configs
from numpy import linalg as la

N0 = 12
S = 1
N_range = np.array([N0])
e_grounds = []      #
eprime_grounds = [] # For sign problems
gap_grounds = []

eplus = 0
eminus = 0
e = 0
for N in N_range:
    print('SIMULATION PARAMS')
    print(N)
    S = int(N0)
    H = sphere.sphere_Hamiltonian(N=N, M=2*S+1, S=S)
    H.generate_basis()
    #H.show_basis()
    H.construct_Hamiltonian()
    #H.print_matrix(H.many_body_H)
    evalues, evecs = H.diagonalise()
    if (N == N0+1):
        eplus = H.e_ground/N
    if (N == N0-1):
        eminus = H.e_ground/N
    if (N == N0):
        e = H.e_ground/N
        e_grounds.append(H.e_ground)
        #print('Hamiltonian eigenvalues ')
        #print(evalues)
        eprime_grounds.append(H.check_sign_problem())
        gap_grounds.append(N0*(eplus + eminus - 2*e))

print('N0 = ', N0)
print('S = ', S)
print(e_grounds)
print(eprime_grounds)
print(gap_grounds)
        
