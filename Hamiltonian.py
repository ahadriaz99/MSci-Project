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
        
        self.many_body_H = np.zeros((self.fock_size(), self.fock_size()))
        
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
        if (fock_basis1.is_vacuum() or fock_basis2.is_vacuum()):
            return 0
        
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
        #print('Basis 1 orderings: ', product1)
        #print('Basis 2 orderings: ', product2)

        normalisation = 1/np.sqrt(product1*product2)
        
        # Form single-body basis overlap matrix and calculate it's permanent
        # Exploit orthonormality of single-particle basis states
        # Optimisation needed for sparse SP_overlaps!
        SP_overlap = np.zeros((self.M, self.M))
        for a in fock_basis1.occup_basis:
            for b in fock_basis2.occup_basis:
                    SP_overlap[a,b] = Dirac_Delta(a, b)
                    
        #print(SP_overlap)
        permanent = ryser.perm_ryser(SP_overlap)
        
        return normalisation*permanent
    
    def generate_basis(self):
        '''
        Generate many-body basis states from repeated combinations
        and index them
        '''
        config_input = np.array(configs.configurations(self.N, self.M)) # Calculate repeated combinations
        assert len(config_input) == self.fock_size() # Check dimensionality

        for i in range(len(config_input)):
            assert len(config_input[i]) == self.M # Check correct input format
            self.basis.append(fock_vector(self.N, self.M, config_input[i],index=i)) # Create fock vectors
            print(config_input[i])

            
    def act_H(self, basis, i, j, k, l):
        '''
        Act Hamiltonian b_i^dagger b_j^daggar b_k b_l on some fock basis state
        '''
        print('Applying Hamiltonian in ', i, j, k, l)
        result_basis = None
        result_basis_1 = None
        total_prefactor = 1
        if i in basis.occup_basis:
            result_basis_1, prefactor_1 = basis.annihilate(i)
            total_prefactor *= prefactor_1
            result_basis, prefactor_2 = result_basis_1.creation(i)
            total_prefactor *= prefactor_2
        else:
            return (fock_vector(self.N, self.M, np.zeros(self.M)), 0)
        
        return (result_basis, total_prefactor)
        
    def matrix_overlap_disc(self, V0, i, j, k, l):
        '''
        Construct many-body matrix elements for disc Hamiltonian
        '''
        return 1#Dirac_Delta(i+j, k+l)#V0*math.factorial(i+l)/2**(i+j)*Dirac_Delta(i+j, k+l)
        
        
    def construct_Hamiltonian(self):
        '''
        Construct many-body Hamiltonian explicitly -- OPTIMISATION!
        '''
        assert len(self.basis) == self.fock_size() # Check if basis generation has been invoked
        for basis1 in self.basis:
            for basis2 in self.basis:
                matrix_overlap = 1
                new_basis2, prefactor = self.act_H(basis1, 0, 0, 0, 0)
                print('Many-body indices: ',basis1.index, basis2.index)
                print('Index 1 Fock vector')
                basis1.print_info()
                print('Index 2 Fock vector')
                basis2.print_info()
                print('Resulting Fock vector')
                new_basis2.print_info()
                print('Prefactor ', prefactor)
                
                print('Overlap of basis1 and new_basis2 ', self.overlap(basis1, new_basis2))
                
                
                self.many_body_H[basis1.index, basis2.index] = \
                        matrix_overlap*self.overlap(basis1, new_basis2)
                '''
                for i in range(self.M):
                    for j in range(self.M):
                        for k in range(self.M):
                            for l in range(self.M):
                                matrix_overlap = self.matrix_overlap_disc(1, i, j, k, l)
                                if (matrix_overlap != 0.0):
                                    # Create basis vector from acting Hamiltonian
                                    new_basis2 = self.act_H(basis2, i, j, k, l)
                                    #assert np.sum(new_basis2[0].occups.values()) == self.N # Assert boson number conservation
                                    # Assign matrix element with matrix_overlap including total prefactors and overlaps
                                    #assert (self.overlap(basis1, new_basis2[0])!= 0.0)
                                    self.many_body_H[basis1.index, basis2.index] += \
                                        matrix_overlap*new_basis2[1]*self.overlap(basis1, new_basis2[0])
                  '''
    def print_Hamiltonian(self):
        print(self.many_body_H)
        return
    
    def diagonalise(self):
        '''
        Carry out exact diagonalisation
        '''
                    
def Disc_BEC(Hamiltonian):
    def __init__(self, N, M):
        super().__init__(N, M)
                

fock_basis1 = fock_vector(5, 3, [3, 1, 1])
fock_basis2 = fock_vector(5, 3, [2, 2, 1])
H = Hamiltonian(5, 3)

print(H.overlap(fock_basis1, fock_basis2))
print(H.overlap(fock_basis2, fock_basis1))
#print(H.fock_size())
#H.generate_basis()
#H.construct_Hamiltonian()
#H.print_Hamiltonian()
