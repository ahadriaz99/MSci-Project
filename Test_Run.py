# -*- coding: utf-8 -*-
""""
Testing script for identifying incompressible liquid phases at half-filling

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

N0 = 3
M = N0*(N0-1)
mu = N0**2/(2*M)
N_range = np.array([N0-1,N0, N0+1])
e_grounds = []      
e_values = []
eprime_grounds = [] # For sign problems
L_range = np.linspace(0,M,M+1)

for L in L_range:
    print('SIMULATION PARAMS')
    print(L)
    H = disc.disc_Hamiltonian(N=N0,M = M, L = L)
    #H = sphere.sphere_Hamiltonian(N=N0, M=2*S+1, S=S, L = L)
    H.generate_basis()
    #H.show_basis()
    H.construct_Hamiltonian()
    #H.print_matrix(H.many_body_H)
    evalues, evecs = H.diagonalise()
    e_values.append(evalues)
    e_grounds.append(H.e_ground)
    eprime_grounds.append(H.check_sign_problem())    

print('N0 = ', N0)
print(e_grounds)
print(eprime_grounds)
#print(e_values)

plt.rcParams.update({'figure.figsize': [12, 12], 'xtick.labelsize' : 20,
'ytick.labelsize' : 20, 'font.size': 20, 'lines.markersize': 20} )

#Plotting the energy spectrum 
plt.figure()
plt.grid()
plt.title('Disc Geometry - Contact Interaction with N = {0} and M = {1}'.format(N0,M))
plt.xlabel('Angular Momentum L')
plt.ylabel('Energy Spectrum')
for i,j in enumerate(e_values):
    #print(j)
    #print([L_range[i]]*len(j))
    plt.plot([L_range[i]]*len(j), j, '_', markersize = 20, mew =3)
plt.show()   

#Ground state energy spectrum 
plt.figure()
plt.grid()
plt.title('N = {0}, M = {1}'.format(N0, M))
plt.plot(L_range,e_grounds,'x', label = 'Hamiltonian', markersize = 12)
plt.plot(L_range, eprime_grounds, 'o', label = 'SP Hamiltonian', markersize = 8)
plt.xlabel('Angular Momentum L')
plt.ylabel('Ground State Energy')
plt.legend()
plt.show()


        
