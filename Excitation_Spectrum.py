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

N0 = 5
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

#print('N0 = ', N0)
#print(e_grounds)
#print(eprime_grounds)
print(e_values)

with open('Disc_Full_Spectrum_N%d_M%d.txt'%(N0, M), 'a') as f:
    f.write('N %d M %d \n'%(N0, M))
    for L in L_range:
        f.write('L %d \n'%(L))
        for i in range(len(e_values[int(L)])):
            f.write(str(e_values[int(L)][i])+" ")
        f.write('\n')

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
plt.figure(1)
plt.title('Disc Geometry (Contact repulsion) \n No. bosons N = {0} No. Landau levels M = {1}'.format(N0,M))
plt.xlabel('Total angular momentum L')
plt.ylabel('Energy spectrum [$V_0$]')
for i,j in enumerate(e_values):
    #print(j)
    #print([L_range[i]]*len(j))
    plt.plot([L_range[i]]*len(j), j, '_', markersize = 25, mew =5)
plt.grid()
plt.savefig('Disc_Full_Spectrum_N%d_M%d.jpeg'%(N0, M))

#Ground state energy spectrum
plt.figure(2)
plt.title('Disc Geometry (Contact repulsion) \n No. bosons N = {0} No. Landau levels M = {1}'.format(N0,M))
plt.plot(L_range,e_grounds,'x', label = 'Hamiltonian', markersize = 15, color='red', mew=3)
plt.plot(L_range, eprime_grounds, 'o', label = 'SP Hamiltonian', markersize = 15, color='blue')
plt.xlabel('Total angular momentum L')
plt.ylabel('Ground state energy [$V_0$]')
plt.legend()
plt.grid()
plt.savefig('Disc_Ground_State_N%d_M%d.jpeg'%(N0, M))

