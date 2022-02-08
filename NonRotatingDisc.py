# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 16:09:50 2022

@author: Owner
"""
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
    
class non_rotatingHamiltonian(Hamiltonian):
    
    '''
    Core class for impelementing Full-Configuration Interaction
    Hamiltonians in bosonic Fock space basis for non-rotating pancakes
    '''
    def __init__(self,N, M, L=0):
        super().__init__(N,M)
        self.tolerance = 1e-10
        self. L = L # Restrict total angular momentum for each Fock vector
        
        self.V0 = 1
    
    def generate_basis(self):
        '''
        Generate many-body basis states from repeated combinations
        and index them
        '''
        config_input = np.array(configs.configurations(self.N, self.M)) # Calculate repeated combinations
        #assert len(config_input) == self.fock_size # Check dimensionality
        
        index = 0
        for i in range(len(config_input)):
            assert len(config_input[i]) == self.M # Check correct input format
            self.basis.append(fock_vector(self.N, self.M, config_input[i],index=i, S=self.S))
            index += 1
        
        print('Basis generation complete')
        print('Fock space size: ', self.fock_size)
                
        self.basis = np.array(self.basis)
        self.fock_size = index
        self.many_body_H = np.zeros((self.fock_size, self.fock_size))
        print(self.fock_size)
    
    
    def matrix_overlap(self, i, j, k, l):
        '''
        Construct many-body overlap matrix for disc Hamiltonian
        '''
        
        self.additionalfactor = 1
        
        if (i+j != k+l):
            return 0
        else:
            return self.V0*self.additionalfactor
    
    def diag_entry(self, basis):
        '''
        Returns diagonal entry for contact repulsion Hamiltonian
        '''
        assert len(self.basis) == self.fock_size # Check if basis generation has not been invoked
        diag_element = 0
        
        occup_basis = np.sort(basis.occup_basis)
        #print(basis.print_info())
        #print(occup_basis)
        for index in range(len(occup_basis)):
            i = occup_basis[index]
            if basis.occups[i] > 1:
                #print(i)
                # Half factor comes from Hamiltonian definition
                diag_element += 0.5*self.matrix_overlap(i, i, i, i)\
                                *basis.occups[i]*(basis.occups[i]-1)
            # we only have to consider non-equal i, j pairs as symmetry
            # gives equal elements for ijij jiij, ijji, jiji basis indices
 
            for jndex in range(index+1, len(occup_basis)):
                j = occup_basis[jndex]
                #print(i, j)

                diag_element += 2*self.matrix_overlap(i, j, i, j)\
                                *basis.occups[i]*(basis.occups[j])
        return diag_element

    def construct_off_diag_entries(self, basis):
        
        off_diag_element = 0
        occup_basis = np.sort(basis.occup_basis)
        new_occups = np.zeros(self.M)
        for index in range(len(occup_basis)):
            i = occup_basis[index]
            
            for jndex in range(index, len(occup_basis)):
                j = occup_basis[jndex]
                for k in range(self.M):
                    new_occups = np.zeros(self.M)
                    new_basis_index = None
                    l = i + j - k
                    if (l >= self.M or l < 0):
                        continue
                    if (k == i or k == j or l == k or l == i or l == j):
                        continue
                    if (i != j):
                        # Copy basis occupation
                        for q in basis.occup_basis:
                            new_occups[q] = basis.occups[q]
                        # Construct basis with non-zero entry
                        new_occups[i] = basis.occups[i] - 1
                        new_occups[j] = basis.occups[j] - 1
                        if (k in basis.occups):
                            new_occups[k] = basis.occups[k] + 1
                        else:
                            new_occups[k] =  1
                            
                        if (l in basis.occups):
                            new_occups[l] = basis.occups[l] + 1
                        else:
                            new_occups[l] = 1
                            
                        new_fock = fock_vector(self.N, self.M, new_occups)
                        new_basis_index = None
                        # Search newly constructed basis index
                        for basis2 in self.basis:
                            if basis2.occups == new_fock.occups:
                                if (basis2.index != basis.index):
                                    new_basis_index = basis2.index
                                    break
                        if (new_basis_index is None):
                            print('New basis not in Hamiltonian space')
                            print(new_fock.print_info)
                            self.show_basis()
                            assert 0
                        
                        # Assign matrix element
                        self.many_body_H[basis.index, new_basis_index] = \
                        2*np.sqrt(basis.occups[i]*basis.occups[j]*new_occups[k]*new_occups[l])*self.matrix_overlap(i, j, k, l)
                        self.many_body_H[new_basis_index, basis.index] = self.many_body_H[basis.index, new_basis_index]
                        
                    else:
                        if (basis.occups[i] < 2):
        
                            continue
                        # Construct basis with non-zero entry for i = j
                        for q in basis.occup_basis:
                            new_occups[q] = basis.occups[q]
                        # See Wilkin paper for angular momentum transfer rules
                        new_occups[i] = basis.occups[i] - 2
                        if (k in basis.occups):
                            new_occups[k] = basis.occups[k] + 1
                        else:
                            new_occups[k] =  1
                            
                        if (l in basis.occups):
                            new_occups[l] = basis.occups[l] + 1
                        else:
                            new_occups[l] = 1
                            
                        new_fock = fock_vector(self.N, self.M, new_occups)
                        new_basis_index = None
                        # Search newly constructed basis index
                        for basis2 in self.basis:
                            if basis2.occups == new_fock.occups:
                                if (basis2.index != basis.index):
                                    new_basis_index = basis2.index
                                    break
                        if (new_basis_index is None):
                            print('New basis not in Hamiltonian space')
                            print(new_fock.print_info)
                            self.show_basis()
                            assert 0
                        
                        # Assign matrix element
                        self.many_body_H[basis.index, new_basis_index] = \
                        np.sqrt(basis.occups[i]*(basis.occups[i]-1)*new_occups[k]*new_occups[l])*self.matrix_overlap(i, j, k, l)
                        self.many_body_H[new_basis_index, basis.index] = self.many_body_H[basis.index, new_basis_index]
                            
        
    def construct_Hamiltonian_fast(self):
        # Wilkin exact eigenstates paper prescription
        assert len(self.basis) == self.fock_size # Check if basis generation has not been invoked
    
        # Diagonal entries
        #print(self.basis)
        print('Hamiltonian construction...')
        print('Fock size: ', self.fock_size)
        counter = 1
        for basis in self.basis:
            self.many_body_H[basis.index, basis.index] = self.diag_entry(basis)
            self.construct_off_diag_entries(basis)
            if (counter % 100 == 0):
                print('Fast Hamiltonian construction progress [%] ', (counter/self.fock_size)*100)
            counter += 1
            
H = non_rotatingHamiltonian(2,2)
H.generate_basis()
H.construct_Hamiltonian_fast()
H.print_matrix(H.many_body_H)

evalues, evecs = H.diagonalise()
print('Hamiltonian eigenvalues [V0]')
print(evalues)
print('Ground state energy [V0] ', H.e_ground)
H.check_sign_problem()


    
    
