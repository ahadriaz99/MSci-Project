# -*- coding: utf-8 -*-
"""
Finite Size Scaling Script 

by Ahad Riaz and Marcell Kovacs
"""

from Hamiltonian import Hamiltonian
import math
import matplotlib.pyplot as plt
import copy
from scipy.optimize import curve_fit
import numpy as np
import sympy as sym

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

def FiniteSize(Length_ratio):
    
    Nstart, Nend = 2, 10
    
    #S0 = 4
    #M0 = 2*S0 + 1
    #Length_ratio = 1e-6
    N = np.linspace(int(Nstart),int(Nend),int((Nend-Nstart+1)))
    N_inverse = []
    for i in range(Nstart, Nend+1, 2):
        N_inverse.append(1/i)
    print(N_inverse)
    print(N)
    
    MF_amps = []
    GP_amps = []
    condensate_fractions = []
    depletion_fractions = []
    
    for i in range(Nstart, Nend+1, 2):
        H = non_rotatingHamiltonian(N=i,S= i/2,M= i + 1, length_ratio = Length_ratio)
        H.generate_basis()
        H.construct_Hamiltonian_fast()
        
        evalues, evecs = H.diagonalise()
        print('Hamiltonian eigenvalues [V0]')
        print(evalues)
        print('Ground state energy [V0] ', H.e_ground)
        print('Ground state configuration', H.e_vector_ground)
        #H.show_basis()
        MF_amp, GP_amp, condensate_fraction = H.ground_state_analysis()
        MF_amps.append(MF_amp)
        GP_amps.append(GP_amp)
        condensate_fractions.append(condensate_fraction)
        depletion_fractions.append(1-condensate_fraction)
        #H.check_sign_problem()
        #H.check_degeneracy()
    
    return condensate_fractions, depletion_fractions

Nstart, Nend = 2, 10
N = np.linspace(int(Nstart),int(Nend),int((Nend-Nstart+1)))
N_inverse = []

for i in range(Nstart, Nend+1, 2):
    N_inverse.append(1/i)
print(N_inverse)


def func(x, a, b, c):
    return a*np.log(b*x) + c

CF1, DF1 = FiniteSize(1e-10)
CF2, DF2 = FiniteSize(1e-1)
CF3, DF3 = FiniteSize(5e-1)
CF4, DF4 = FiniteSize(1e0)
CF5, DF5 = FiniteSize(1.5e0)
CF6, DF6 = FiniteSize(2e0)

popt, pcov = curve_fit(func, N_inverse, CF5, p0 = [0.05,0.005,1])
popt1, pcov1 = curve_fit(func, N_inverse, CF2, p0 = [0.05,0.005,1])
print ("a = %s , b = %s, c = %s" % (popt[0], popt[1], popt[2]))
#plt.plot(N_inverse, func(N_inverse, *popt), label="Fitted Curve") #same as line above \/

plt.figure()
plt.title('Finite Size Scaling for Condensate Fraction \n' +\
          'N = %d to %d '%(Nstart, Nend))
plt.ylabel('Condensate Fraction (Log Scale)')
plt.xlabel('1/N (Log Scale)')
plt.grid()
plt.plot(np.log(N_inverse), np.log(CF1), color = 'black', label = 'LR = 1e-10')
plt.plot(np.log(N_inverse), np.log(CF1), 'x', color = 'black', mew = 3)
plt.plot(np.log(N_inverse), np.log(CF2), color = 'red', label = 'LR = 1e-1')
plt.plot(np.log(N_inverse), np.log(CF2), 'x', color = 'red', mew = 3)
plt.plot(np.log(N_inverse), np.log(CF3), color = 'blue', label = 'LR = 5e-1')
plt.plot(np.log(N_inverse), np.log(CF3), 'x', color = 'blue', mew = 3)
plt.plot(np.log(N_inverse), np.log(CF4), color = 'orange', label = 'LR = 1e0')
plt.plot(np.log(N_inverse), np.log(CF4), 'x', color = 'orange', mew = 3)
plt.plot(np.log(N_inverse), np.log(CF5), color = 'green', label = 'LR = 1.5e0')
plt.plot(np.log(N_inverse), np.log(CF5), 'x', color = 'green', mew = 3)
plt.plot(np.log(N_inverse), np.log(CF6), color = 'purple', label = 'LR = 2e0')
plt.plot(np.log(N_inverse), np.log(CF6), 'x', color = 'purple', mew = 3)
#plt.plot(np.linspace(0,0.5,100), popt[0]*np.log(popt[1]*np.linspace(0,0.5,100)) + popt[2], label="Fitted Curve")
plt.legend()
plt.show()

plt.figure()
plt.title('Finite Size Scaling for Condensate Fraction \n' +\
          'N = %d to %d '%(Nstart, Nend))
plt.ylabel('Condensate Fraction')
plt.xlabel('1/N')
plt.grid()
plt.plot(N_inverse, CF1, color = 'black', label = 'LR = 1e-10')
plt.plot(N_inverse, CF1, 'x', color = 'black', mew = 3)
plt.plot(N_inverse, CF2, color = 'red', label = 'LR = 1e-1')
plt.plot(N_inverse, CF2, 'x', color = 'red', mew = 3)
plt.plot(N_inverse, CF3, color = 'blue', label = 'LR = 5e-1')
plt.plot(N_inverse, CF3, 'x', color = 'blue', mew = 3)
plt.plot(N_inverse, CF4, color = 'orange', label = 'LR = 1e0')
plt.plot(N_inverse, CF4, 'x', color = 'orange', mew = 3)
plt.plot(N_inverse, CF5, color = 'green', label = 'LR = 1.5e0')
plt.plot(N_inverse, CF5, 'x', color = 'green', mew = 3)
plt.plot(N_inverse, CF6, color = 'purple', label = 'LR = 2e0')
plt.plot(N_inverse, CF6, 'x', color = 'purple', mew = 3)
plt.legend()
plt.show()

plt.figure()
plt.title('Finite Size Scaling for Condensate Fraction \n' +\
          'N = %d to %d and Length Ratio = 1.5'%(Nstart, Nend))
plt.ylabel('Condensate Fraction')
plt.xlabel('1/N')
plt.grid()
plt.plot(N_inverse, CF5)
plt.plot(N_inverse, CF5, 'x', color = 'black', mew = 3)
#plt.plot(np.linspace(0,0.5,100), popt[0]*np.log(popt[1]*np.linspace(0,0.5,100)) + popt[2], label="Fitted Curve")
plt.legend()
plt.show()

plt.figure()
plt.title('Finite Size Scaling for Condensate Fraction \n' +\
          'N = %d to %d and Length Ratio = 0.1'%(Nstart, Nend))
plt.ylabel('Condensate Fraction')
plt.xlabel('1/N')
plt.grid()
plt.plot(N_inverse, CF2)
plt.plot(N_inverse, CF2, 'x', color = 'black', mew = 3)
#plt.plot(np.linspace(0,0.5,100), popt[0]*np.log(popt[1]*np.linspace(0,0.5,100)) + popt[2], label="Fitted Curve")
plt.legend()
plt.show()
        