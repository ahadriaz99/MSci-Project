# -*- coding: utf-8 -*-
""""
SubClass to generate the many body Hamiltonian for sphere geometry

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
    
class sphere_Hamiltonian(Hamiltonian):
    
    def __init__(self,N,M,S):
        '''Additional argument for angular momentum S'''
        super().__init__(N,M)
        self.S = S
        self.tolerance = 1e-10
          
    def matrix_overlap_sphere(self, i, j, k, l):
        '''
        Construct many-body matrix elements for disc Hamiltonian
        '''
        S = self.S
        sum_SM = 0
        
        sum_SM += math.factorial(S + i) * math.factorial(S-i)
        sum_SM += math.factorial(S + j) * math.factorial(S-j)
        sum_SM += math.factorial(S + k) * math.factorial(S-k)
        sum_SM += math.factorial(S + l) * math.factorial(S-l)
        
        print(sum_SM)
        
        
        V0 = 1
        return Dirac_Delta(i+j, k+l)*V0*((math.factorial(2*S+1))**2 * math.factorial(2*S + i + j) *
                                         math.factorial(2*S - i - j))/(S * math.factorial(4*S + 1) * np.sqrt(sum_SM))
        
    def H_element(self, basis1, basis2):
        '''
        Calculate matrix element between 2 many-body basis states
        '''
        #assert basis1 in self.basis and basis2 in self.basis
        #print('Basis 1')
        #basis1.print_info()
        #print('Basis 2')
        #basis2.print_info()
        element = 0 # Matrix element
        # Loop over all possible single-particle state indices
        # NEEDS OPTIMISATION
        for i in range(self.M):
            for j in range(self.M):
                for k in range(self.M):
                    for l in range(self.M):
                    
                        #print('Operator indices', i, j, k, l)
                        #print('BRA BASIS')
                        #basis1.print_info()
                        #print('KET BASIS')
                        #basis2.print_info()
                        
                        # Calculate scalar integral overlaps
                        matrix_overlap = self.matrix_overlap_sphere(i, j, k, l)
                        if (matrix_overlap != 0):
                            # Apply annihlation on basis 1 BRA -> (i, j) must be reversed
                            new_basis1, total_prefactor_1 = self.annihilate(basis1, j, i)
                            # Apply annihilation on basis 2 KET
                            new_basis2, total_prefactor_2 = self.annihilate(basis2, k, l)
                            #print('TRANSFORMED BASIS1')
                            #new_basis1.print_info()
                            #print('TRANSFORMED BASIS2')
                            #new_basis2.print_info()
                            #print('Prefactor ', total_prefactor_1*total_prefactor_2)
                            # Assemble matrix element contribution
                            element += 0.5*matrix_overlap*self.overlap(new_basis1, new_basis2)*total_prefactor_1*total_prefactor_2
                            #print('Overlap: ', self.overlap(new_basis1, new_basis2))
        return element

H = sphere_Hamiltonian(5,3,4)

H.generate_basis()
H.show_basis()
H.construct_Hamiltonian()
H.print_matrix(H.many_body_H)
