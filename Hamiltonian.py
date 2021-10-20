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
        Overlap of many-body basis states 1 and 2 given that SP states are
        orthonormal
        '''
        
        if (fock_basis1.is_vacuum() or fock_basis2.is_vacuum()):
            return 0
        
        # First calculate orderings of excitations for
        # each basis state
        
        print('Basis1', fock_basis1.occup_basis.sort())
        print('Basis2', fock_basis2.occup_basis.sort())

        product = 1
        # NOT WORKING
        if (fock_basis1.occups.keys() != fock_basis2.occups.keys()):
            return 0
        else:
        # We have ensured occupied bases are the same
            for i in fock_basis1.occup_basis:
                # For orthonormal SP states, basis mismatch renders overlap to zero
                # Calculate orderings of bosons in SP basis states
                product *= math.factorial(fock_basis1.occups[i])*\
                           math.factorial(fock_basis2.occups[i])
            
            return 1/np.sqrt(product)
        

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
            print(config_input[i], i)
            
    def basis_overlap(self)
        '''
        Form overlap between basis states --> See if basis is non-orthogonal,
        if the returne matrix is non-diagonal
        '''
        
        assert self.basis != []
        
        basis_overlap = np.array((self.fock_size(), self.fock_size()))
        
        for basis1 in self.basis:
            for basis2 in self.basis:
                basis_overlap[basis1.index, basis2.index] = self.overlap(basis1, basis2)
        
        return basis_overlap
        
    def basis_orthogonalisation(self):
        '''
        Create new basis vectors by some orthogonalisation procedure
        '''
            
    def act_H(self, basis, i, j, k, l):
        '''
        Act Hamiltonian b_i^dagger b_j^daggar b_k b_l on some fock basis state
        '''
        #print('Applying Hamiltonian in ', i, j, k, l)
        result_basis = None
        result_basis_1 = None
        result_basis_2 = None
        result_basis_3 = None
        prefactor_1 = 1
        prefactor_2 = 1
        prefactor_3 = 1
        prefactor_4 = 1


        total_prefactor = 1
        '''
        # ONLY FOR TESTING
        if i in basis.occup_basis:
            result_basis_1, prefactor_1 = basis.annihilate(i)
            total_prefactor *= prefactor_1
            result_basis, prefactor_2 = result_basis_1.creation(i)
            total_prefactor *= prefactor_2
        else:
            return (fock_vector(self.N, self.M, np.zeros(self.M)), 0)
        '''
        assert np.array([i, j, k, l]).all() < self.M and np.array([i, j, k, l]).all() >= 0

        vacuum = fock_vector(self.N, self.M, np.zeros(self.M))
        if ((i+j) != (k+l)):
            return (vacuum, 0)
        if (l not in basis.occup_basis):
            return (vacuum, 0)
        result_basis_1, prefactor_1 = basis.annihilate(l)
        if (k not in result_basis_1.occup_basis):
            return (vacuum, 0)
        result_basis_2, prefactor_2 = result_basis_1.annihilate(k)
        result_basis_3, prefactor_3 = result_basis_2.creation(j)
        result_basis, prefactor_4 = result_basis_3.creation(i)
        total_prefactor = prefactor_1*prefactor_2*prefactor_3*prefactor_4
        
        #print(sum(result_basis.occups.values()), self.N)
        assert sum(result_basis.occups.values()) == self.N
        
        
        return (result_basis, total_prefactor)
        
    def matrix_overlap_disc(self, i, j, k, l):
        '''
        Construct many-body matrix elements for disc Hamiltonian
        '''
        V0 = 1
        return Dirac_Delta(i+j, k+l)#*V0*math.factorial(i+l)/2**(i+j)
        
    def H_element(self, basis1, basis2):
        assert basis1 in self.basis and basis2 in self.basis
        print('Basis 1')
        basis1.print_info()
        print('Basis 2')
        basis2.print_info()
        element = 0
        for i in range(self.M):
            for j in range(self.M):
                for k in range(self.M):
                    for l in range(self.M):
                        
                        print('Apply operators ', i, j, k, l)
                        print('BRA BASIS')
                        basis1.print_info()
                        
                        print('KET BASIS')
                        basis2.print_info()
                        
                        matrix_overlap = self.matrix_overlap_disc(i, j, k, l)
                    
                        new_basis2, total_prefactor = self.act_H(basis2, i, j, k, l)
                        print('TRANSFORMED BASIS')
                        new_basis2.print_info()
                        print('Prefactor ', total_prefactor)
                        element += matrix_overlap*self.overlap(basis1, new_basis2)*total_prefactor
                        print('Overlap: ', self.overlap(basis1, new_basis2))
        return element
        
    def construct_Hamiltonian(self):
        '''
        Construct many-body Hamiltonian explicitly -- OPTIMISATION!
        '''
        assert len(self.basis) == self.fock_size() # Check if basis generation has been invoked
        for basis1 in self.basis:
            for basis2 in self.basis:
                self.many_body_H[basis1.index, basis2.index] = self.H_element(basis1, basis2)

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
                
N = 2
M = 3
H = Hamiltonian(N, M)


                

#print(H.fock_size())
H.generate_basis()
H.construct_Hamiltonian()
H.print_Hamiltonian()
#H.H_element(H.basis[0], H.basis[2])

