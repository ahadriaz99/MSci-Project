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

params = {
   'axes.labelsize': 35,
   'font.size': 35,
   'legend.fontsize': 35,
   'lines.linewidth' : 2,
   'lines.markersize' : 35,
   'xtick.labelsize': 35,
   'ytick.labelsize': 35,
   'figure.figsize': [40, 20]
   }
plt.rcParams.update(params)

N = 3
S = 2
L_range = np.linspace(0, N*S, N*S+1)
e_grounds = []
eprime_grounds = []

for L in L_range:
    print('SIMULATION PARAMS')
    print('N=',N, 'S=', S, 'L=',L)
    H = sphere.sphere_Hamiltonian(N=N, M=2*S+1, S=S, L=L)
    H.generate_basis()
    #H.show_basis()
    H.construct_Hamiltonian()
    #H.print_matrix(H.many_body_H)
    evalues, evecs = H.diagonalise()
    e_grounds.append(evalues.min())
    eprime_grounds.append(H.check_sign_problem())
    
plt.title('N = %d, S = %d, L = %d'%(N, S, N*S))
plt.plot(L_range, e_grounds, 'x', label='Hamiltonian')
plt.plot(L_range, eprime_grounds, '.', label='SP Hamiltonian')
plt.xlabel('Total angular momentum L')
plt.ylabel('Ground state energy')
plt.grid()
plt.legend()
plt.savefig('Excitations_N%dS%d.jpeg'%(N, S))

