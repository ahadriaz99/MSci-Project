'''
Class to generate bosonic Hamiltonians in occupation number representation

14/10/21

by Marcell Dorian Kovacs and Ahad Riaz
'''

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
        
class Hamiltonian:
    '''
    Core class for impelementing Full-Configuration Interaction
    Hamiltonians in bosonic Fock space basis
    '''
    def __init__(self, N, M, S=0):
    
        self.tolerance = 1e-8 # Numerical tolerance on equalities...etc.
        
        self.N = N # Total number of bosons/excitations
        self.M = M # Total number of single particle basis states
        self.S = S # Added for spherical Hamiltonian
        
        # Size of Hilbert space (no. repeated combinations available)
        self.fock_size =  int(math.factorial(self.N+self.M-1)/\
               (math.factorial(self.N)*math.factorial(self.M-1)))
        
        # Dictionary indexed by many-body basis state index for state A, B,
        # storing the overlap of the matrix
        # Many-body basis state index will come from generating combinations
        self.overlaps = {}
        
        # Array of Fock basis vectors forming Fock space
        self.basis = []
        
        # Many-body Hamiltonian -> OPTIMISATION
        self.many_body_H = np.zeros((self.fock_size, self.fock_size))
        
        # Ground state properties
        self.e_ground = None
        self.e_vector_ground = None
        
    
    def overlap(self, fock_basis1, fock_basis2):
        '''
        Overlap of many-body basis states 1 and 2 given that SP states are
        orthonormal
        '''
        
        if fock_basis1.is_vacuum() and fock_basis2.is_vacuum():
            return 1
        elif fock_basis1.is_vacuum() or fock_basis2.is_vacuum():
            return 0
        
        if fock_basis1.occups != fock_basis2.occups:
            return 0
        else:
            return 1
            #print(fock_basis1.occups)
            #for basis in fock_basis1.occups.values():
            #    #print(basis)
            #    overlap *= 1/math.factorial(basis)
        #return overlap
                

    def generate_basis(self):
        '''
        Generate many-body basis states from repeated combinations
        and index them
        '''
        config_input = np.array(configs.configurations(self.N, self.M)) # Calculate repeated combinations
        assert len(config_input) == self.fock_size # Check dimensionality

        for i in range(len(config_input)):
            assert len(config_input[i]) == self.M # Check correct input format
            self.basis.append(fock_vector(self.N, self.M, config_input[i],index=i, S=self.S)) # Create fock vectors
            #print(config_input[i], i)
            
    def basis_overlap(self):
        assert (self.basis != [])
        basis_overlap = np.zeros((self.fock_size, self.fock_size))
        for basis1 in self.basis:
            for basis2 in self.basis:
                basis_overlap[basis1.index, basis2.index] = self.overlap(basis1, basis2)
        return basis_overlap
    
    def show_basis(self):
        print('Many-body configurations')
        print('------------------------')
        for basis in self.basis:
            basis.print_info()
        print('------------------------')
        
    def annihilate(self, basis, i, j):
        '''
        Apply annihilation on a KET vector on indices i, j
        '''
        assert np.array([i, j]).all() < self.M and np.array([i, j]).all() >= 0
        
        result_basis_1 = None
        result_basis = None
        
        prefactor_1 = 1
        prefactor_2 = 1

        vacuum = fock_vector(self.N, self.M, np.zeros(self.M))
        if (j not in basis.occup_basis):
            return (vacuum, 0)
        result_basis_1, prefactor_1 = basis.annihilate(j)
        if (i not in result_basis_1.occup_basis):
            return (vacuum, 0)
        result_basis, prefactor_2 = result_basis_1.annihilate(i)
        
        return (result_basis, prefactor_1*prefactor_2)
        
    #def matrix_overlap_disc(self, i, j, k, l):
        #'''
        #Construct many-body matrix elements for disc Hamiltonian
        #'''
        #V0 = 1
        #return Dirac_Delta(i+j, k+l)*V0*math.factorial(i+j)/2**(i+j)/\
        #       np.sqrt(math.factorial(i)*math.factorial(j)*math.factorial(k)*math.factorial(l))
        
    def H_element(self, basis1, basis2):
        '''
        Calculate matrix element between 2 many-body basis states
        '''
        #assert basis1 in self.basis and basis2 in self.basis
        #print('Basis 1')
        #basis1.print_info()
        #print('Basis 2')
        #basis2.print_info()
        #element = 0 # Matrix element
        # Loop over all possible single-particle state indices
        # NEEDS OPTIMISATION
        #for i in range(self.M):
        #    for j in range(self.M):
        '''
            for k in range(self.M):
                    for l in range(self.M):
                    
                        #print('Operator indices', i, j, k, l)
                        #print('BRA BASIS')
                        #basis1.print_info()
                        #print('KET BASIS')
                        #basis2.print_info()
                        
                        # Calculate scalar integral overlaps
                        matrix_overlap = self.matrix_overlap_disc(i, j, k, l)
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
        return element'''
        
    def construct_Hamiltonian(self):
        '''
        Construct many-body Hamiltonian explicitly -- OPTIMISATION NEEDED!
        '''
        
        #print(len(self.basis), (self.fock_size))
        assert len(self.basis) == self.fock_size # Check if basis generation has been invoked
        i = 0
        for basis1_index in range(len(self.basis)):
            # Exploit hermitianity
            for basis2_index in range(basis1_index, len(self.basis)):
                element = self.H_element(self.basis[basis1_index], self.basis[basis2_index])
                self.many_body_H[basis1_index, basis2_index] = element
                self.many_body_H[basis2_index, basis1_index] = element
                    
                i += 2
                if (i % 100 == 0):
                    print('Hamiltonian construction progress [%]', (i/(self.fock_size**2))*100)
        
        return self.many_body_H
                
    def print_matrix(self, matrix):
        
        isDiagonal = True
        isSymmetric = True
    
        print('Many-body matrix')
        print('---------------------')
        for basis1 in self.basis:
            for basis2 in self.basis:
                print("%2.4f"%matrix[basis1.index, basis2.index], end="         ")
                if (abs(matrix[basis1.index, basis2.index] - matrix[basis2.index, basis1.index]) > self.tolerance):
                    isSymmetric = False
                if abs(matrix[basis1.index, basis2.index]) > self.tolerance and basis1.index != basis2.index:
                    isDiagonal = False
            print()
        print('Diagonal? ', isDiagonal)
        print('Symmetric? ', isSymmetric)
        print('---------------------')
        assert (isSymmetric)

    def diagonalise(self):
        '''
        Carry out exact diagonalisation to produce ground state energy
        and configuration
        '''
        evalues, evectors = la.eigh(self.many_body_H)
        e_ground = evalues.min()
        
        self.e_ground, self.e_vector_ground = e_ground, evectors[np.where(evalues == e_ground)]
        
        return evalues, evectors
        
    def check_sign_problem(self):
        '''
        Diagonalise H' = -(1-delta_{ij})|H_{ij}| + H_{ij}delta_{ij}
        transformed Hamiltonian to check severity of Monte Carlo sign problem
        Compare lowest eigenvalues of H and H' -- for moderate sign problem,
        the difference should be small
        '''
        assert self.e_ground is not None
        H_prime = -abs(self.many_body_H)
        H_prime[np.diag_indices(self.fock_size)] = np.diag(self.many_body_H)
        #self.print_matrix(H_prime)
        
        evalues, evectors = la.eigh(H_prime)
        e_ground_prime = evalues.min()
        
        print('Sign problem analysis')
        print('---------------------')
        print(self.e_ground, 'Ground state of phsyical H [V0]')
        print(e_ground_prime, 'Ground state of unphysical H\' [V0]')
        
        return e_ground_prime
        
        
        
def test_overlap():
    N = 5
    M = 4
    H = Hamiltonian(N, M)
    vector1 = fock_vector(N, M, [4, 1, 0, 0])
    vector2 = fock_vector(N, M, [2, 3, 0, 0])
    vector3 = fock_vector(N, M, [2, 0, 1, 2])
    vector4 = fock_vector(N, M, [4, 1, 0, 0])
    print(H.overlap(vector1, vector2))
    print(H.overlap(vector1, vector3))
    print(H.overlap(vector1, vector4))
    print(H.overlap(vector2, vector3))
                
#N = 3
#M = 6
#H = Hamiltonian(N, M)

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
