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
    
class sphere_Hamiltonian_fast(Hamiltonian):
    
    def __init__(self,N,M,S,L=-1):
        '''Additional argument for angular momentum S'''
        assert (2*S+1 == M)
        self.S = S
        self.M = 2*self.S + 1

        super().__init__(N,M,S)

        self.tolerance = 1e-10
        self.L = L
        self.V0 = 1
        #self.v = np.zeros((self.M), self.M):
        #for index in range(self.M):
         #   for jndex in range(self.M):
         ##       i = index - self.S
          #      j = jndex - self.S
          #      self.v[index, jndex] = np.sqrt(float(math.factorial(self.S+(i+j))*math.factorial(self.S-(i+j))))/\
           #                    np.sqrt(float(math.factorial(self.S + i)*math.factorial(self.S - i)*\
           #                                  math.factorial(self.S - j)*math.factorial(self.S - j)))*\
           #                    np.sqrt(self.S*math.factorial(4*self.S + 1))*float(math.factorial(2*self.S+1))
        
          
    def matrix_overlap_sphere(self, i, j, k, l):
        '''
        Construct many-body matrix elements for sphere Hamiltonian
        i, j, k, l refer to values between [-S, S]
        '''
        if (i+j != k+l):
            return 0
        else:
            S = self.S
            sum_SM = 1
        
            sum_SM *= math.factorial(S + i) * math.factorial(S - i)
            sum_SM *= math.factorial(S + j) * math.factorial(S - j)
            sum_SM *= math.factorial(S + k) * math.factorial(S - k)
            sum_SM *= math.factorial(S + l) * math.factorial(S - l)
        
            #print(sum_SM)
            
            return (self.V0*((math.factorial(2*S+1))**2 * math.factorial(2*S + i + j) *
                            math.factorial(2*S - i - j)))/(S * math.factorial(4*S + 1) * np.sqrt(float(sum_SM)))
            #return self.V0*self.v[i + self.S, j + self.S]*self.v[k + self.S, l + self.S]
        
    def generate_basis(self):
        '''
        Generate many-body basis states from repeated combinations
        and index them
        '''
        # If basis generation has not been invoked, generate basis
        print('Basis generation...')
        configs.sphere_config_fast(int(self.N), int(self.M), int(self.L), int(self.S))
        index = 0
        file='Sphere_Configurations_N%dM%dL%dS%d.txt'%(self.N, self.M, self.L, self.S)
        print('Reading in configurations...')
        with open(file, 'r') as f:
            for line in f:
                split_line = line.split()
                if (split_line[0]=='N,'):
                    continue
                N = split_line[0]
                M = split_line[1]
                L = split_line[2]
                S = split_line[3]
                basis = []
                config = split_line[4:]
                #print(config)
                for item in config:
                    #print(item)
                    basis.append(int(item))
                #print(N, M, L, len(config), len(basis))
                #print(config, basis)
                #print(self.N, N, self.M, M)
                assert int(N) == self.N and int(M) == self.M and int(M) == 2*(self.S) +1
                vector = fock_vector(int(N), int(M), np.array(basis), int(S))
                assert vector.ang_mom() == self.L
                vector = fock_vector(int(N), int(M), np.array(basis), S= int(S), index=index)
                self.basis.append(vector)
                index += 1
                if (index % 100 == 0):
                    print('No. basis read-in ', index)
            
        print('Basis generation complete')
        print('Fock space size: ', self.fock_size)
                
        self.basis = np.array(self.basis)
        self.fock_size = index
        self.many_body_H = np.zeros((self.fock_size, self.fock_size))
        #print(config_input[i], i)
          
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
            #print('occups basis', basis.occups[i])
            if basis.occups[i] > 1:
                #print(i)
                # Half factor comes from Hamiltonian definition
                #print('i',i - self.S)
                #print(self.matrix_overlap_sphere(i - self.S, i- self.S, i- self.S, i- self.S))
                diag_element += 0.5*self.matrix_overlap_sphere(i- self.S, i- self.S, i- self.S, i- self.S)\
                                *basis.occups[i]*(basis.occups[i]-1)
                #print(diag_element)
            # we only have to consider non-equal i, j pairs as symmetry
            # gives equal elements for ijij jiij, ijji, jiji basis indices
 
            for jndex in range(index+1, len(occup_basis)):
                j = occup_basis[jndex]
                #print(i, j)
                #print('i,j',i - self.S, j- self.S)
                #print(self.matrix_overlap_sphere(i - self.S, j- self.S, i- self.S, j- self.S))
                #print(basis.occups[i])
                #print(basis.occups[j])
                diag_element += 2*self.matrix_overlap_sphere(i - self.S, j- self.S, i- self.S, j- self.S)\
                                *basis.occups[i]*(basis.occups[j])
                #print(diag_element)
        return diag_element
    
    def construct_off_diag_entries(self, basis):
            
        off_diag_element = 0
        occup_basis = np.sort(basis.occup_basis)
        new_occups = np.zeros(self.M)
        for index in range(len(occup_basis)):
            i = occup_basis[index]
            #print('here is i', i)
            for jndex in range(index, len(occup_basis)):
                j = occup_basis[jndex]
                #print('here is j', j)
                for k in range(self.M):
                    new_occups = np.zeros(self.M)
                    new_basis_index = None
                    l = i + j - k
                    if (l >= self.M-1 or l < 0):
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
                        2*np.sqrt(basis.occups[i]*basis.occups[j]*new_occups[k]*new_occups[l])*self.matrix_overlap_sphere(i - self.S, j- self.S, k- self.S, l- self.S)
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
                        np.sqrt(basis.occups[i]*(basis.occups[i]-1)*new_occups[k]*new_occups[l])*self.matrix_overlap_sphere(i- self.S, j- self.S, k- self.S, l- self.S)
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
                        if (i+j != k+l):
                            continue
                    
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
#
#H = sphere_Hamiltonian_fast(N=2, M=3, S=1, L=0)
#H.generate_basis()
#print('Fock size', H.fock_size)
#H.construct_Hamiltonian_fast()
##H.show_basis()
#H.print_matrix(H.many_body_H)
