# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 14:09:25 2021

@author: Owner
"""

from Hamiltonian import Hamiltonian

import numpy as np
import math


#H = ham.Hamiltonian(2,3)


#H.generate_basis()
#H.show_basis()
#H.construct_Hamiltonian()
#H.print_matrix(H.many_body_H)

#H.show_basis()
#H.print_matrix(H.basis_overlap())

#evalues, evecs = H.diagonalise()
#print('Hamiltonian eigenvalues [V0]')
#print(evalues)
#print('Ground state energy [V0] ', H.e_ground)

#H.check_sign_problem()

def Dirac_Delta(a, b):
    '''
    Simple delta function
    '''
    if (a == b):
        return 1
    else:
        return 0
    
class sub_Hamiltonian(Hamiltonian):
    
    def __init__(self,N,M):
        super().__init__(N,M)
        self.tolerance = 1e-10
          
    def matrix_overlap_disc(self, i, j, k, l):
        '''
        Construct many-body matrix elements for disc Hamiltonian
        '''
        V0 = 1
        return Dirac_Delta(i+j, k+l)*V0*math.factorial(i+j)/2**(i+j)/\
               np.sqrt(math.factorial(i)*math.factorial(j)*math.factorial(k)*math.factorial(l))
        
    def H_element(self, basis1, basis2):
        '''
        Calculate matrix element between 2 many-body basis states
        '''
        assert basis1 in self.basis and basis2 in self.basis
        #print('Basis 1')
        #basis1.print_info()
        #print('Basis 2')
        #basis2.print_info()
        element = 0 # Matrix element
        # Loop over all possible single-particle state indices
        # NEEDS OPTIMISATION

        return element

class BoseHubbard1D(Hamiltonian):
    def __init__(self,N,M, t, U, mu):
        super().__init__(N,M)
        
        self.t = t
        self.U = U
        self.mu = mu
        
    def H_element(self, basis1, basis2):
        element = 0
        for i in range(self.M):
            if (i < self.M-1):
                j = i+1
            else:
                j = 0
            
            new_basis1, pref1 = basis1.annihilate(i)
            new_basis2, pref2 = basis2.annihilate(j)
            
            element += -self.t*self.overlap(new_basis1, new_basis2)*pref1*pref2
            
            new_basis1, pref1 = basis1.annihilate(j)
            new_basis2, pref2 = basis2.annihilate(i)
            
            element += -self.t*self.overlap(new_basis1, new_basis2)*pref1*pref2
            
            if (i in basis2.occup_basis):
                element += -self.mu*basis2.occups[i]*self.overlap(basis1, basis2)
                element += self.U/2*(basis2.occups[i])*(basis2.occups[i]-1)*self.overlap(basis1, basis2)
        
        return element

    def one_particle_spectrum(self):
            '''
            Construct Bloch waves as a non-interacting spectrum
            '''
            n = np.linspace(0, self.M-1, self.M)
            print('n ', n)
            H_analytic = -2*self.t*np.array(np.cos(2*np.pi*n/self.M))
            return H_analytic
        
H = BoseHubbard1D(2, 5, t=1, U=0, mu=0)


H.generate_basis()
H.show_basis()
H.construct_Hamiltonian()
H.print_matrix(H.many_body_H)
evalues, evecs = H.diagonalise()
print('Hamiltonian eigenvalues ')
print(evalues)
print(evecs[np.where(evalues.min())])
print(H.one_particle_spectrum())

#H.print_matrix(H.basis_overlap())
#H.show_basis()
#H.print_matrix(H.basis_overlap())
