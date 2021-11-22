# -*- coding: utf-8 -*-
"""
Excitation Spectrum for Spherical Geometry

25/10/21

by Marcell Dorian Kovacs and Ahad Riaz
"""

from Hamiltonian import Hamiltonian
import numpy as np
import math
import matplotlib.pyplot as plt
import copy
import timeit

from Fock_vector import fock_vector
from Sphere_Hamiltonian import sphere_Hamiltonian
from Sphere import sphere_Hamiltonian_fast
import Ryser_Algorithm as ryser
import config as configs
from numpy import linalg as la

N0 = 10
S = 4
flux = 2*S
M = 2*S + 1
mu = N0/(2*S)
N_range = np.array([N0-1,N0, N0+1])
e_grounds = []
e_values = []
eprime_grounds = [] # For sign problems
gaps = []
L_range = np.linspace(0,M,M+1)

print('Basis generation...')
configs.sphere_config_fast(N0, M, M, S)
print('Basis generation complete!')

for L in L_range:
    print('SIMULATION PARAMS')
    print('L ', L)
    H = sphere_Hamiltonian_fast(N=N0, M = M, L = L, S=S)
    H.generate_basis()
    #H.show_basis()
    print('Fock space ', H.fock_size)
    if (H.fock_size > 5000):
        print('Large Fock space...')
        #assert 0
        

    H.construct_Hamiltonian_fast()
    #H.print_matrix(H.many_body_H)
    evalues, evecs = H.diagonalise()
    e_values.append(evalues)
    e_grounds.append(H.e_ground)
    eprime_grounds.append(H.check_sign_problem())
    gaps.append(H.gap)
    
    with open('Sphere_Full_Spectrum_N%d_M%d_S%d.txt'%(N0, M, S), 'a') as f:
        if (L == 0):
            f.write('N %d M %d \n'%(N0, M))
        f.write('L %d \n'%(L))
        for i in range(len(e_values[int(L)])):
            f.write(str(e_values[int(L)][i])+" ")
        f.write('\n')
            
#print('N0 = ', N0)
#print(e_grounds)
#print(eprime_grounds)

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
plt.title('Sphere Geometry (Contact repulsion) \n No. bosons N = {0} No. Landau levels M = {1} Flux Quanta 2S = {2}'.format(N0,M,flux))
plt.xlabel('Total angular momentum L')
plt.ylabel('Energy spectrum [$V_0$]')
for i,j in enumerate(e_values):
    #print(j)
    #print([L_range[i]]*len(j))
    plt.plot([L_range[i]]*len(j), j, '_', markersize = 25, mew =5)
plt.grid()
plt.savefig('Sphere_Full_Spectrum_N%d_M%d_S%d.jpeg'%(N0, M,S))

#Ground state energy spectrum
plt.figure(2)
plt.title('Sphere Geometry (Contact repulsion) \n No. bosons N = {0} No. Landau levels M = {1} Flux Quanta 2S = {2}'.format(N0,M,flux))
plt.plot(L_range,e_grounds,'x', label = 'Hamiltonian', markersize = 15, color='red', mew=3)
plt.plot(L_range, eprime_grounds, 'o', label = 'SP Hamiltonian', markersize = 15, color='blue')
plt.xlabel('Total angular momentum L')
plt.ylabel('Ground state energy [$V_0$]')
plt.legend()
plt.grid()
plt.savefig('Sphere_Ground_State_N%d_M%d_S%d.jpeg'%(N0, M, S))
plt.figure(3)
plt.title('Sphere Geometry (Contact repulsion) \n No. bosons N = {0} No. Landau levels M = {1} Flux Quanta 2S = {2}'.format(N0,M,flux))
plt.plot(L_range,gaps,'x', label = 'Hamiltonian', markersize = 15, color='red', mew=3)
plt.ylabel('Local gap $|E_0 - E_1|$ [$V_0$]')
plt.legend()
plt.grid()
plt.savefig('Sphere_Local_Gap_N%d_M%d_S%d.jpeg'%(N0, M, S))