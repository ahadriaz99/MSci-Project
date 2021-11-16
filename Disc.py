# -*- coding: utf-8 -*-
"""
SubClass to generate the many body Hamiltonian for disc geometry

25/10/21

by Marcell Dorian Kovacs and Ahad Riaz
"""
from Hamiltonian import Hamiltonian
import numpy as np
import math
import matplotlib.pyplot as plt
import copy

from Fock_vector import fock_vector
from Disc_Hamiltonian import disc_Hamiltonian
import Ryser_Algorithm as ryser
import config as configs
from numpy import linalg as la
    
class disc_Hamiltonian_fast(Hamiltonian):
    
    def __init__(self,N, M, L):
        super().__init__(N,M)
        self.tolerance = 1e-10
        self. L = L # Restrict total angular momentum for each Fock vector
        
        self.V0 = 1
        self.v = np.zeros((self.M, self.M))
        for i in range(self.M):
            for j in range(i, self.M):
                #print(i, j)
                ifac = math.factorial(i)
                jfac = math.factorial(j)
                self.v[i, j] = np.sqrt(float(math.factorial(i+j)))/\
                               np.sqrt(float(ifac*jfac))
                self.v[j, i] = self.v[i, j]
                
    def large_fac(self, I, J):
        max = np.max([I, J])
        min = np.min([I, J])
        result = 1
        divisor = max
        for j in range(abs(I-J), I+J+1):
            result /= divisor
            if (divisor >= 1):
                divisor -= 1
            result *= j
        if (divisor >= 1):
            result /= math.factorial(divisor)
        return float(result)
        
    def generate_basis(self):
        '''
        Generate many-body basis states from repeated combinations
        and index them
        '''
        # If basis generation has not been invoked, generate basis
        print('Basis generation...')
        configs.disc_config(int(self.N), int(self.M), int(self.L))
        index = 0
        file='Disc_Configurations_N%dM%dL%d.txt'%(self.N, self.M, self.L)
        print('Reading in configurations...')
        with open(file, 'r') as f:
            for line in f:
                split_line = line.split()
                if (split_line[0]=='N,'):
                    continue
                N = split_line[0]
                M = split_line[1]
                L = split_line[2]
                basis = []
                config = split_line[3:]
                for item in config:
                    basis.append(int(item))
                #print(N, M, L, len(config), len(basis))
                #print(config, basis)
                #print(self.N, N, self.M, M)
                assert int(N) == self.N and int(M) == self.M
                vector = fock_vector(int(N), int(M), np.array(basis))
                assert vector.ang_mom() == self.L
                vector = fock_vector(int(N), int(M), np.array(basis), S=0, index=index)
                self.basis.append(vector)
                index += 1
                if (index % 100 == 0):
                    print('Index ', index)
            
        print('Basis generation complete')
        print('Fock space size: ', self.fock_size)
                
        self.basis = np.array(self.basis)
        self.fock_size = index
        self.many_body_H = np.zeros((self.fock_size, self.fock_size))
        #print(config_input[i], i)
          
    def matrix_overlap_disc(self, i, j, k, l):
        '''
        Construct many-body overlap matrix for disc Hamiltonian
        '''
        if (i+j != k+l):
            return 0
        else:
            return self.V0*self.v[i,j]*self.v[k, l]/2**(i+j)
            
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
            #print(i)
            # Half factor comes from Hamiltonian definition
            diag_element += 0.5*self.matrix_overlap_disc(i, i, i, i)\
                            *basis.occups[i]*(basis.occups[i]-1)
            # we only have to consider non-equal i, j pairs as symmetry
            # gives equal elements for ijij jiij, ijji, jiji basis indices
 
            for jndex in range(index+1, len(occup_basis)):
                j = occup_basis[jndex]
                #print(i, j)

                diag_element += 2*self.matrix_overlap_disc(i, j, i, j)\
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
                        2*np.sqrt(basis.occups[i]*basis.occups[j]*new_occups[k]*new_occups[l])*self.matrix_overlap_disc(i, j, k, l)
                        self.many_body_H[new_basis_index, basis.index] = self.many_body_H[basis.index, new_basis_index]
                        
                    else:
                        if (basis.occups[i] < 2):
                            continue
                        # Construct basis with non-zero entry for i = j
                        for q in basis.occup_basis:
                            new_occups[q] = basis.occups[q]
                        
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
                        np.sqrt(basis.occups[i]*(basis.occups[i]-1)*new_occups[k]*new_occups[l])*self.matrix_overlap_disc(i, j, k, l)
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
            if (counter % 10 == 0):
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
                        #print('Operator indices', i, j, k, l)
                        #print('BRA BASIS')
                        #basis1.print_info()
                        #print('KET BASIS')
                        #basis2.print_info()
                        
                        # Calculate scalar integral overlaps
                        matrix_overlap = self.matrix_overlap_disc(i, j, k, l)
                        if (matrix_overlap == 0):
                            continue
                        else:
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
        
    def BoseDump(self):
        '''
        Output matrix elements in the Bose dump format
        '''
        with open('Disc_Overlaps_N%dM%dL%d.txt'%(self.N, self.M, self.L), 'a') as f:
            f.write('&FCI NMODE= %d, NBOSON= %d\n&END\n'%(self.M, self.N))
            for i in range(self.M):
                for j in range(self.M):
                    for k in range(self.M):
                        for l in range(self.M):
                            matrix_overlap = self.matrix_overlap_disc(i, j, k, l)
                            if (abs(matrix_overlap) > self.tolerance):
                                f.write(('%5.10f %d %d %d %d \n'%(matrix_overlap, (i+1), (j+1), (k+1), (l+1))))
        f.close()
                                
#N = 5
#M = 5
#L = 5
#H = disc_Hamiltonian_fast(N=N,M=M,L=L)
#H.generate_basis()
#H.construct_Hamiltonian_fast()
#H.print_matrix(H.many_body_H)
#H1 = disc_Hamiltonian(N=N,M=M,L=L)
#H1.generate_basis()
#H.show_basis()
#H1.construct_Hamiltonian()

#H.show_basis()
#H1.show_basis()

#H.print_matrix(H.many_body_H)
#H1.print_matrix(H1.many_body_H)
#H.print_matrix(H.many_body_H)

#configs.configurations(N=10, M=10)
#H.generate_basis()
#H.show_basis()
#H.BoseDump()
#H.print_matrix(H.construct_Hamiltonian())
#evalues, evecs = H.diagonalise()
#print('Hamiltonian eigenvalues [V0]')
##print(evalues)
#print('Ground state energy [V0] ', H.e_ground)
#H.check_sign_problem()


#print(configs.configurations(N=3, M=10))
#H.generate_basis()
#H.show_basis()
#H.construct_Hamiltonian()
#H.print_matrix(H.many_body_H)
#evalues, evecs = H.diagonalise()#
#print('Hamiltonian eigenvalues [V0]')
#print(evalues)
#print('Ground state energy [V0] ', H.e_ground)
#H.check_sign_problem()
#>>>>>>> 5913d15 (Minor updates)

