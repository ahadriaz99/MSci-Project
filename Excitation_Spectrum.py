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

N = 6
M = 3
filling = N**2/(2*M)
L_max = 10
L_range = np.linspace(0, L_max, L_max+1)
e_grounds = []
eprime_grounds = []
evalues_all = []

for L in L_range:
    print('SIMULATION PARAMS')
    print('N=',N, 'M=', M, 'L=',L)
    H = disc.disc_Hamiltonian(N=N, M=M, L=L)
    H.generate_basis()
    #H.show_basis()
    H.construct_Hamiltonian()
    #H.print_matrix(H.many_body_H)
    evalues, evecs = H.diagonalise()
    e_grounds.append(evalues.min())
    evalues_all.append(evalues)
    eprime_grounds.append(H.check_sign_problem())
    
plt.title('Disc Hamiltonian N = %d, M = %d, L = %d'%(N, M, N*(M-1)))
plt.plot(L_range, e_grounds, 'x', label='Hamiltonian')
plt.plot(L_range, eprime_grounds, '.', label='SP Hamiltonian')
plt.xlabel('Total angular momentum L')
plt.ylabel('Ground state energy')
plt.grid()
plt.legend()
plt.savefig('Disc_Excitations_N%dM%d.jpeg'%(N, M))

print('Filling')
print(filling)
print('Eigenvalues')
for L_index in range(len(L_range)):
    print('L= ', L_range[L_index])
    print('Evalues')
    print(evalues_all[L_index])
