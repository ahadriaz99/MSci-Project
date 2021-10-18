'''
Class to generate bosonic Hamiltonians in occupation number representation

14/10/21

by Marcell Dorian Kovacs
'''

import numpy as np
import math
import matplotlib.pyplot as plt
import copy

from Fock_vector import fock_vector
import Ryser_Algorithm as ryser
import config as configs

def Dirac_Delta(a, b):
    '''
    Simple delta function
    '''
    if (a == b):
        return 1
    else:
        return 0
        
class Hamiltonian:
    '''
    Core class for impelementing Full-Configuration Interaction
    Hamiltonians in bosonic Fock space basis
    '''
    def __init__(self, N, M):
        
        self.N = N # Total number of bosons/excitations
        self.M = M # Total number of single particle basis states
        
        # Dictionary indexed by many-body basis state index for state A, B,
        # storing the overlap of the matrix
        # Many-body basis state index will come from generating combinations
        self.overlaps = {}
        
        # Array of Fock basis vectors forming Fock space
        self.basis = []
        
    
    def fock_size(self):
        '''
        Number of many-body basis states formed by repeated combinations
        Choose N excitations repeatadly on M basis states
        '''
        return int(math.factorial(self.N+self.M-1)/\
               (math.factorial(self.N)*math.factorial(self.M-1)))
        
    
    def overlap(self, fock_basis1, fock_basis2):
        '''
        Overlap of many-body basis states 1 and 2
        '''
        
        # First calculate orderings of excitations for
        # each basis state
        
        # Orderings of basis1 occupations
        product1 = 1
        for item in fock_basis1.occups.values():
            product1 *= math.factorial(item)
        # Orderdings of basis2 occupations
        product2 = 1
        for item in fock_basis2.occups.values():
            product2 *= math.factorial(item)
        print('Basis 1 orderings: ', product1)
        print('Basis 2 orderings: ', product2)

        normalisation = 1/np.sqrt(product1*product2)
        
        # Form single-body basis overlap matrix and calculate it's permanent
        # Exploit orthonormality of single-particle basis states
        # Optimisation needed for sparse SP_overlaps!
        SP_overlap = np.zeros((self.M, self.M))
        for a in fock_basis1.occup_basis:
            for b in fock_basis2.occup_basis:
                    SP_overlap[a,b] = Dirac_Delta(a, b)
                    
        print(SP_overlap)
        permanent = ryser.perm_ryser(SP_overlap)
        
        return normalisation*permanent
    
    def generate_basis(self):
        '''
        Generate many-body basis states from repeated combinations
        and index them
        '''
        config_input = np.array(configs.configurations(self.N, self.M))
        assert len(config_input) == self.fock_size()

        for config in config_input:
            self.basis.append(fock_vector(self.N, self.M, config))
        for basis in self.basis:
            basis.print_info()
        
    def generate_overlaps(self):
        '''
        Construct many-body matrix elements specifically for Hamiltonians
        implemented in subclasses
        '''
        
        
def Disc_BEC(Hamiltonian):
    def __init__(self, N, M):
        super().__init__(N, M)
                

#fock_basis1 = fock_vector(5, 3, [3, 1, 1])
#fock_basis2 = fock_vector(5, 3, [3, 1, 1])
H = Hamiltonian(2, 3)
H.generate_basis()
