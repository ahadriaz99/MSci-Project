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
from Sphere import sphere_Hamiltonian_fast
import config as configs
from numpy import linalg as la




N0 = 7
S = N0-1
M = 2*S + 1
L = S
mu = N0/(2*S)


#N_range = np.array([N0-1,N0, N0+1])
e_grounds = []      
e_values = []
eprime_grounds = [] # For sign problems
L_range = np.linspace(0,L,L+1)
gap = []

eplus = 0
eminus = 0 
e = 0

for N in range(N0+1):
    print('SIM PARAMS')
    print('N ', N)
    N_range = np.array([N-1,N, N+1])
    if N == 0 or N == 1 or N==2 or N==3:
        gap.append(0)
    else:
        for i in N_range:
            print('Gap loop')
            print('N, M, L, S', i, M, L, S)
            H = sphere_Hamiltonian_fast(i, M, S, L)
            H.generate_basis()
            H.construct_Hamiltonian_fast()
            evalues,evecs = H.diagonalise()
            if i == N + 1:
                eplus = H.e_ground/i
            if i == N - 1:
                eminus = H.e_ground/i
            if i == N:
                e = H.e_ground/i
                                        
        gap.append(N*(eplus + eminus - 2*e))

with open('Sphere_Gap_Spectrum_Laughlin.txt', 'a') as f:
    for i in range(N0+1):
        f.write('N '+str(i))
        f.write('\n')
        f.write(str(gap[i])+'\n')
        
#Plotting the energy spectrum
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

plt.figure(1)
plt.title('Sphere Geometry (Contact repulsion) \n No. Landau levels M = 2*S + 1  \n Half filling')
plt.ylabel('Energy gap [$V_0$]')
plt.xlabel('$N_0$')
plt.plot(range(N0+1), gap, 'x', color='red')
plt.grid()
plt.savefig('Sphere_Gap_Spectrum_Laughlin1_N%d.jpeg'%(N0))

plt.figure(2)
plt.title('Sphere Geometry (Contact repulsion) \n No. Landau levels M = 2*S + 1  \n Half filling')
plt.ylabel('Energy gap [$V_0$]')
plt.xlabel('1/$N_0$')
plt.plot(1/np.array(range(N0+1))[1:], gap[1:], 'x', color='red')
plt.grid()
plt.savefig('Sphere_Gap_Spectrum_Laughlin2_N%d.jpeg'%(N0))
