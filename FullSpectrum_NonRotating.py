# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 16:44:11 2022

@author: Owner
"""
from Hamiltonian import Hamiltonian
import numpy as np
import math
import matplotlib.pyplot as plt
import copy
import NonRotatingDisc as rd

from Fock_vector import fock_vector
import Ryser_Algorithm as ryser
import Sphere_Hamiltonian as sphere
import Disc as disc
import config as configs
from numpy import linalg as la

N0 = 5
M0 = 10
#mu = N0**2/(2*M)
N_range = np.array([N0-1,N0, N0+1])
e_grounds = []
e_values = []
eprime_grounds = [] # For sign problems
gaps = []
#L_range = np.linspace(0,M-1,M)
M_range = np.linspace(1,M0,M0)
print(M_range)

for M in M_range:
    print('SIMULATION PARAMS')
    print('M ', M)
    H = rd.non_rotatingHamiltonian(N=int(N0), M = int(M))
    #H = sphere.sphere_Hamiltonian(N=N0, M=2*S+1, S=S, L = L)
    H.generate_basis()
    #H.show_basis()
    print('Fock space ', H.fock_size)
    if (H.fock_size > 10000):
        print('Large Fock space...')
        assert 0
    H.construct_Hamiltonian_fast()
    #H.print_matrix(H.many_body_H)
    evalues, evecs = H.diagonalise()
    e_values.append(evalues)
    e_grounds.append(H.e_ground)
    eprime_grounds.append(H.check_sign_problem())
    gaps.append(H.gap)

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

#Plotting the energy spectrum
plt.figure()
plt.title('Disc Geometry (Contact repulsion) \n No. bosons N = {0} No. Landau levels M = {1}'.format(N0,M))
plt.xlabel('Total angular momentum L')
plt.ylabel('Energy spectrum [$V_0$]')
for i,j in enumerate(e_values):
    #print(j)
    #print([L_range[i]]*len(j))
    plt.plot([M_range[i]]*len(j), j, '_', markersize = 25, mew =5)
plt.grid()

#Ground state energy spectrum
plt.figure()
plt.title('Disc Geometry (Contact repulsion) \n No. bosons N = {0} No. Landau levels M = {1}'.format(N0,M))
plt.plot(M_range,e_grounds,'x', label = 'Hamiltonian', markersize = 15, color='red', mew=3)
#plt.plot(M_range, eprime_grounds, 'o', label = 'SP Hamiltonian', markersize = 15, color='blue')
plt.xlabel('Total angular momentum L')
plt.ylabel('Ground state energy [$V_0$]')
plt.legend()
plt.grid()


