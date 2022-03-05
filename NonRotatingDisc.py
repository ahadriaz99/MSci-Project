# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 16:09:50 2022

@author: Owner
"""
from Hamiltonian import Hamiltonian
import numpy as np
import math
import matplotlib.pyplot as plt
import copy

from Hamiltonian import Hamiltonian
from Fock_vector import fock_vector
import Ryser_Algorithm as ryser
import config as configs
from numpy import linalg as la
from decimal import Decimal
from Fock_vector import fock_vector
from Disc_Hamiltonian import disc_Hamiltonian

params = {
   'axes.labelsize': 35,
   'font.size': 35,
   'legend.fontsize': 35,
   'lines.linewidth' : 2,
   'lines.markersize' : 35,
   'xtick.labelsize': 35,
   'ytick.labelsize': 35,
   'figure.figsize': [40, 20]
   }
plt.rcParams.update(params)

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
    def __init__(self,N, M, S, length_ratio=10, L=0):
        super().__init__(N,M)
        self.tolerance = 1e-10
        self. L = L # Restrict total angular momentum for each Fock vector
        self.S = S
        assert M == 2*S + 1
        self.M = M
        # Set interaction energy scale to 1
        self.V0 = 1
        self.lengthratio = length_ratio # = (a_z/a_s: trap length/scattering length)
        # Scale kinetic energy scale accordingly
        self.T0 = np.sqrt(np.pi/2)*np.pi*self.lengthratio
        
        self.condensate_fraction = None # No. excitations in lowest SP state
        self.GP_amplitude = None # Weight of Gross-Pitaevskii (fully condensed) permanent
        self.GP_index = None
        self.MF_perm = None # Amplitude of dominant permanent in FCI expansion
        self.MF_energy = None # Energy content of dominant permanent in FCI expansion
        self.MF_amplitude = None # Amplitude of Gross-Pitaevskii permanent
        
    def generate_basis(self):
        '''
        Generate many-body basis states from repeated combinations
        and index them
        '''
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
                if (self.S in vector.occup_basis):
                    if (vector.occups[self.S] == self.N):
                        self.GP_index = index
                index += 1
                if (index % 100 == 0):
                    print('No. basis read-in ', index)
            
        print('Basis generation complete')
        print('Fock space size: ', self.fock_size)
                
        self.basis = np.array(self.basis)
        self.fock_size = index
        self.many_body_H = np.zeros((self.fock_size, self.fock_size))
    
    
    def matrix_overlap(self, i, j, k, l):
        '''
        Construct many-body overlap matrix for disc Hamiltonian
        '''
        
        self.additionalfactor = 1
        
        if (i+j != k+l):
            return 0
        else:
            return self.V0*self.additionalfactor
    
    def kineticterm(self, i):
        
        return self.T0*((i-self.S)**2)
    
    
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
            diag_element += self.kineticterm(i)
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
            
    def ground_state_analysis(self):
        # Index of MF permanent
        #print(np.max(self.e_vector_ground.T))
        #print(self.e_vector_ground[0])
        max_index = np.where(max(self.e_vector_ground[0], key=abs) == self.e_vector_ground[0])[0][0]
        print('max index', max_index)
        
        
        self.MF_perm = self.basis[max_index]
        self.MF_amplitude = self.e_vector_ground[0][max_index]
        print('Mean-field permanent info')
        print(self.MF_perm.print_info())
        print('Amplitude: ', self.MF_amplitude)
        
        self.GP_amplitude = self.e_vector_ground[0][self.GP_index]
        #print('All permanents:', self.e_vector_ground[0])
        print('Gross-Pitaevskii permanent amplitude: ', self.GP_amplitude)
        
        #print('Permanent energy: ', self.many_body_H[max_index, max_index])
        #print('MF energy / E0: ',  self.many_body_H[max_index, max_index]/self.e_ground)
        
        # Calculating condensate fraction

        self.condensate_fraction = 0

        for index, amplitude in enumerate(self.e_vector_ground.squeeze()):

            if (self.S not in self.basis[index].occup_basis):
                continue
            else:
                self.condensate_fraction += \
                    abs(amplitude)**2 * self.basis[index].occups[self.S]/self.N

        print('Expected condensate fraction: ', self.condensate_fraction)
        print('Condensate depletion: ', 1-self.condensate_fraction)
        
        return self.MF_amplitude, self.GP_amplitude
 
    def check_degeneracy(self):
        '''
        Find degeneracies within spectrum
        '''
        self.degen_evalues = []
        self.degen_evectors = []
        
        for i in range(len(self.evalues)):
            self.degen_evalues.append([i])
            self.degen_evectors.append([self.evectors.T[i]])
            for j in range(i+1, len(self.evalues)):
                if abs(self.evalues[i] - self.evalues[j]) <= self.tolerance:
                    self.degen_evalues[-1].append(j)
                    self.degen_evectors[-1].append(self.evectors.T[j])
        
        degeneracy = np.zeros(len(self.evalues))
        
        for i in range(len(self.evalues)):
            degeneracy[i] = len(self.degen_evalues[i])
            
        plt.title('Degeneracy of spectrum\n'+\
                 'Disc geometry\nN = %d   M = %d   L = %d'%(self.N, self.M, self.L))
        plt.bar(x=np.arange(1, self.fock_size+1), height=degeneracy)
        plt.xlabel('Sorted eigenstate index')
        plt.ylabel('Degeneracy')
        plt.grid()
        plt.legend()
        plt.savefig('BEC_Degeneracy_N%d_M%d_S%d_L%d.jpeg'%(self.N, self.M, self.S, self.L))
        plt.close()
        
        plt.title('Eigenvalue spectrum\n'+\
                 'Disc geometry\nN = %d   M = %d   L = %d'%(self.N, self.M, self.L))
        nT, binsT, patchesT = plt.hist(x=self.evalues, bins=15, color='red',
                            alpha=0.7, rwidth=0.85, label='FCI Spectrum')
        plt.xlabel('Eigenvalues [$V_0$]')
        plt.ylabel('Degeneracy')
        plt.legend()
        plt.grid()
        plt.savefig('BEC_Spectrum_N%d_M%d_S%d_L%d.jpeg'%(self.N, self.M, self.S, self.L))
        plt.close()
        
        assert (self.evalues.min() == self.evalues[0])
        assert (self.fock_size == len(self.evalues))
        
        
        
        for ground_index in range(len(self.degen_evalues[0])):
            print(len(self.degen_evectors[0]), len(self.degen_evalues[0]))
            assert len(self.degen_evectors[0]) == len(self.degen_evalues[0])
            #print(self.degen_evectors[0][ground_index])
            print(len(self.degen_evectors[0][ground_index]))
            print(ground_index)
            plt.figure(ground_index)
            plt.title('Degenerate ground state configuration index %d \n'%(ground_index)+\
                     'Disc geometry\nN = %d   M = %d   L = %d'%(self.N, self.M, self.L))
            plt.bar(x=np.arange(1, self.fock_size+1), height=self.degen_evectors[0][ground_index][0])
            #nT, binsT, patchesT =\
            #plt.hist(x=self.degen_evectors[0][ground_index],bins=self.fock_size, color='red',alpha=0.7, rwidth=0.85, label='Full Configuration')
            plt.xlabel('Many-body basis index')
            plt.ylabel('Amplitude')
            plt.grid()
            plt.legend()
            plt.savefig('BEC_Ground_Config_i%d_N%d_M%d_S%d_L%d.jpeg'%(ground_index, self.N, self.M,self.S, self.L))
            plt.close()
        
                
        #print('Degenerate evector indices')
        #print(self.degen_evalues)
        #print(self.degen_evectors)


H = non_rotatingHamiltonian(N=3,S=1,M=3)
H.generate_basis()
H.construct_Hamiltonian_fast()
#H.print_matrix(H.many_body_H)

evalues, evecs = H.diagonalise()
print('Hamiltonian eigenvalues [V0]')
print(evalues)
print('Ground state energy [V0] ', H.e_ground)
print('Ground state configuration', H.e_vector_ground)
#H.show_basis()
MF_amp, GP_amp = H.ground_state_analysis()
print(MF_amp, GP_amp)
#H.check_sign_problem()
H.check_degeneracy()

    
