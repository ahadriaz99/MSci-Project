# -*- coding: utf-8 -*-
""""
SubClass to generate the many body Hamiltonian for Bose Hubbard model

25/10/21

by Marcell Dorian Kovacs and Ahad Riaz
"""

from Hamiltonian import Hamiltonian

import numpy as np
import math

def Dirac_Delta(a, b):
    '''
    Simple delta function
    '''
    if (a == b):
        return 1
    else:
        return 0
    
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
            
            new_basis2, pref2 = basis2.annihilate(i)
            new_basis1, pref1 = basis1.annihilate(j)
            
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
        
H = BoseHubbard1D(3, 10, t=1, U=0, mu=0)

H.generate_basis()
H.show_basis()
H.construct_Hamiltonian()
H.print_matrix(H.many_body_H)
#evalues, evecs = H.diagonalise()
#print('Hamiltonian eigenvalues ')
#print(evalues)
#print(H.one_particle_spectrum())
#H.print_matrix(H.basis_overlap())
#H.show_basis()
#H.print_matrix(H.basis_overlap())