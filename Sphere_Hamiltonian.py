# -*- coding: utf-8 -*-
""""
SubClass to generate the many body Hamiltonian for sphere geometry

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
import config as configs
from numpy import linalg as la


def Dirac_Delta(a, b):
    '''
    Simple delta function
    '''
    if (a == b):
        return 1
    else:
        return 0
    
class sphere_Hamiltonian(Hamiltonian):
    
    def __init__(self,N,M,S,L=-1):
        '''Additional argument for angular momentum S'''
        assert (2*S+1 == M)
        self.S = S
        self.M = 2*self.S + 1
        

        super().__init__(N,M,S)

        self.tolerance = 1e-10
        self.L = L
          
    def matrix_overlap_sphere(self, i, j, k, l):
        '''
        Construct many-body matrix elements for disc Hamiltonian
        '''
        S = self.S
        sum_SM = 1
        
        sum_SM *= math.factorial(S + i) * math.factorial(S-i)
        sum_SM *= math.factorial(S + j) * math.factorial(S-j)
        sum_SM *= math.factorial(S + k) * math.factorial(S-k)
        sum_SM *= math.factorial(S + l) * math.factorial(S-l)
        
        #print(sum_SM)
        
        
        V0 = 1
        return Dirac_Delta(i+j, k+l)*V0*((math.factorial(2*S+1))**2 * math.factorial(2*S + i + j) *
                                         math.factorial(2*S - i - j))/(S * math.factorial(4*S + 1) * np.sqrt(float(sum_SM)))
        
    def generate_basis(self):
        '''
        Generate many-body basis states from repeated combinations
        and index them
        '''
        config_input = np.array(configs.configurations(self.N, self.M)) # Calculate repeated combinations
        assert len(config_input) == self.fock_size # Check dimensionality
        
        index = 0
        for i in range(len(config_input)):
            #print(i)
            assert len(config_input[i]) == self.M # Check correct input format
            vector = fock_vector(self.N, self.M, config_input[i], S=self.S)
            # Only add restricted ang. mom. bases to the Fock spaces
            if(vector.ang_mom() == self.L):
                self.basis.append(fock_vector(self.N, self.M, config_input[i],index=index, S=self.S)) # Create fock vectors
                #self.basis[-1].print_info()
                index += 1
        self.fock_size = index
        self.many_body_H = np.zeros((self.fock_size, self.fock_size))
        
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
                        # the indices run from 0 to M == 2*S + 1
                        # for obtaining the correct angular momentum (i, j, k, l) must be between -S to S
                        # Convention we take for Fock vector is to
                        matrix_overlap = self.matrix_overlap_sphere(i-self.S, j-self.S, k-self.S, l-self.S)
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
