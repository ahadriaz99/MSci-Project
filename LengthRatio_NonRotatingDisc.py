
# -*- coding: utf-8 -*-
"""
Varying T0 and V0 Script 

by Ahad Riaz and Marcell Kovacs
"""

from Hamiltonian import Hamiltonian
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import sympy as sym
import copy

from Hamiltonian import Hamiltonian
from NonRotatingDisc import non_rotatingHamiltonian
from Fock_vector import fock_vector
import Ryser_Algorithm as ryser
import config as configs
from numpy import linalg as la
from decimal import Decimal
from Fock_vector import fock_vector
from Disc_Hamiltonian import disc_Hamiltonian

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

def LengthRatioScanner(N0):
    #N0 = 2
    S0 = N0/2
    M0 = 2*S0 + 1
    #length_ratio = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]
    length_ratio = np.linspace(1e-10,1,31)
    MF_amps = []
    GP_amps = []
    condensate_fraction = []

    for i in length_ratio:
        print('Length Ratio', i)
        H = non_rotatingHamiltonian(N=N0,S=S0,M=int(M0), length_ratio = i)
        H.generate_basis()
        H.construct_Hamiltonian_fast()
        
        evalues, evecs = H.diagonalise()
        print('Hamiltonian eigenvalues [V0]')
        print(evalues)
        print('Ground state energy [V0] ', H.e_ground)
        print('Ground state configuration', H.e_vector_ground)
        #H.show_basis()
        MF_amp, GP_amp, cond_frac = H.ground_state_analysis()
        MF_amps.append(MF_amp)
        GP_amps.append(GP_amp)
        condensate_fraction.append(cond_frac)
        #H.check_sign_problem()
        #H.check_degeneracy()
    
    return MF_amps, GP_amps, condensate_fraction

length_ratio = np.linspace(1e-10,1,31)
#length_ratio = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]
MF_4, GP_4, CF_4 = LengthRatioScanner(4)
MF_6, GP_6, CF_6 = LengthRatioScanner(6)
MF_8, GP_8, CF_8 = LengthRatioScanner(8)
MF_10, GP_10, CF_10 = LengthRatioScanner(10)
MF_2, GP_2, CF_2 = LengthRatioScanner(2)


Nstart, Nend = 2, 10
N = np.linspace(int(Nstart),int(Nend),int((Nend-Nstart+1)))
N_inverse = []

for i in range(Nstart, Nend+1, 2):
    N_inverse.append(1/i)
print(N_inverse)


def func(x, a,  b, c):
    return  a - b*np.exp(-c*x)

popt, pcov = curve_fit(func, length_ratio, CF_6)

#%%
plt.figure()
plt.title('Comparison of Kinetic Energy to Interaction Strength using Length Ratio \n' +\
          'N = %d to %d, M = 2*S + 1, S = N/2'%(2, 10))
plt.ylabel('Condensate Fraction')
plt.xlabel('Length Ratio')
#plt.xticks(np.arange(0,1,step=31))
plt.grid()
plt.plot(length_ratio, CF_2, label = 'N=2', color = 'blue')
plt.plot(length_ratio, CF_2, 'x', mew = 3, color = 'blue')
plt.plot(length_ratio, CF_4, label = 'N=4', color = 'green')
plt.plot(length_ratio, CF_4, 'x', mew = 3, color = 'green')
plt.plot(length_ratio, CF_6, label = 'N=6', color = 'purple')
plt.plot(length_ratio, CF_6, 'x', mew = 3, color = 'purple')
plt.plot(length_ratio, CF_8, label = 'N=8', color = 'orange')
plt.plot(length_ratio, CF_8, 'x', mew = 3, color = 'orange')
plt.plot(length_ratio, CF_10, label = 'N =10', color = 'red')
plt.plot(length_ratio, CF_10, 'x', mew = 3, color = 'red')
plt.legend()
plt.show()

     
plt.figure()
plt.title('Comparison of Kinetic Energy to Interaction Strength using Length Ratio \n' +\
          'N = 10, M = 2*S + 1, S = N/2')
plt.ylabel('Condensate Fraction')
plt.xlabel('Length Ratio')
#plt.xticks(np.arange(0,1,step=31))
plt.grid()
plt.plot(length_ratio, CF_10, label = 'N=10', color = 'blue')
plt.plot(length_ratio, CF_10, 'x', mew = 3, color = 'blue')
#plt.plot(np.linspace(0,1,100), popt[0] - popt[1]*np.exp(-popt[2]*np.linspace(0,1,100)), label="Fitted Curve")
plt.legend()
plt.show()