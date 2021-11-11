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



N0=9
#N_range = np.array([N0-1,N0, N0+1])
e_grounds = []      
e_values = []
eprime_grounds = [] # For sign problems

gap = []

eplus = 0
eminus = 0 
e = 0

'''for N in range(N0):
    if N == 0 or N == 1:
        continue
    else:'''


for N in range(N0+1):
    print('SIM PARAMS')
    print('N ', N)
    N_range = np.array([N-1,N, N+1])
    if N == 0 or N == 1 or N==2:
        gap.append(0)
    else:
        for i in N_range:
            print('Gap loop')
            print('N, M, L', i, N*(N), N*(N))
            H = disc.disc_Hamiltonian(N = i, M = N, L = N)
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

with open('Disc_Gap_Spectrum_Full.txt', 'a') as f:
    for i in range(N0+1):
        f.write('N '+str(i))
        f.write('\n')
        f.write(str(gap[i])+'\n')
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
#Plotting the gap spectrum
plt.figure(1)
plt.title('Disc Geometry (Contact repulsion) \n No. Landau levels M = N Total ang. mom. L = N')
plt.ylabel('Energy gap [$V_0$]')
plt.xlabel('$N_0$')
plt.plot(range(N0+1), gap, 'x', color='red')
plt.grid()
plt.savefig('Disc_Gap_Spectrum_High_Filling_N%d.jpeg'%N0)
