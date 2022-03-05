# -*- coding: utf-8 -*-
"""
Finite Size Scaling Script 

by Ahad Riaz and Marcell Kovacs
"""

from Hamiltonian import Hamiltonian
import numpy as np
import math
import matplotlib.pyplot as plt
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

Nstart, Nend = 2, 10

S0 = 4
M0 = 2*S0 + 1
Length_ratio = 10
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
    H = non_rotatingHamiltonian(N=i,S= i/2,M= i + 1)
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
    H.check_degeneracy()

plt.figure()
plt.title('Finite Size Scaling for Condensate Fraction \n' +\
          'N = %d to %d , Length Ratio = %d'%(Nstart, Nend, Length_ratio))
plt.ylabel('Condensate Fraction')
plt.xlabel('1/N')
plt.grid()
plt.plot(N_inverse, condensate_fractions)
plt.plot(N_inverse, condensate_fractions, 'x', color = 'black', mew = 3)
plt.show()

plt.figure()
plt.title('Finite Size Scaling for Depletion Fraction \n' +\
          'N = %d to %d , Length Ratio = %d'%(Nstart, Nend, Length_ratio))
plt.ylabel('Depletion Fraction')
plt.xlabel('1/N')
plt.grid()
plt.plot(N_inverse, depletion_fractions)
plt.plot(N_inverse, depletion_fractions, 'x', color = 'black', mew = 3)
plt.show()
        