'''
Class to generate bosonic Hamiltonians in occupation number representation

13/10/21

by Marcell Dorian Kovacs and Ahad Riaz
'''

import numpy as np
import math
import matplotlib.pyplot as plt
import copy

def Dirac_Delta(a, b):
    '''
    Simple delta function
    '''
    print(a, b)
    if (a == b):
        return 1
    else:
        return 0

class fock_vector:
    def __init__(self, N, M, occupations, index=None):
        '''
        Class for representing a bosonic Fock space vectors in
        occupation number representation
        N: number of particles in Fock space
        M: number of single-particle basis states
        occupations: M size array giving explicit occupation for
            each SP basis
        '''

        self.N = N
        self.M = M
        # Vector representing the occupations of each SP state
        # as dictionary
        self.occups = {}
        # Vector representing which basis states are occupied
        # as an array. Each element represents a basis state index [1,M]
        # which is occupied in the vector
        self.occup_basis = []
        
        self.index = index # Unique index of basis state assigned when creating object
        
        # Check that input occupancy ket is indeed an element of Fock space
        
        # Valid occupation numbers
        assert np.all(np.array(occupations)) >= 0
        #assert type(np.all(np.array(occupations))) == 'int'
        # Correct space of basis states
        assert len(np.array(occupations)) == M
        
        if (np.array(occupations).sum() == 0):
            return
        else:
            # Boson number conservation
            assert np.array(occupations).sum() == N
        
        #Optimised dictionary construction
                
        occups = []
        
        for i in range(len(occupations)):
            # Only store non-zero occupations
            if (occupations[i] != 0):
                # Append occupied basis index
                self.occup_basis.append(i)
                # Store index and occupation
                occups.append((i, occupations[i]))
                
        # Create dictionary of occupations
        # indexed by occupied single particle indices
        
        # Memory saving for long Fock vectors
        self.occups = dict(np.array(occups))

    def is_vacuum(self):
        return self.occup_basis == []
    def annihilate(self, index):
        '''
        Method to carry out annihilation on a Fock ket
        index: basis state on which to act annihilation
        Method does NOT modify original fock vector
        Returns:
            result_vec: Fock vector
            prefactor: annihilation factor
        '''
        assert (index >= 0 and index < self.M)
                
        if (index not in self.occup_basis):
            # No excitations on sought basis
            # Return vacuum
            return (fock_vector(self.N, self.M, np.zeros(self.M)), 0)
        
        result_vec = copy.deepcopy(self)
        # Removing last boson from basis
        if (self.occups[index] == 1):
            # Remove occupation item
            del result_vec.occups[index]
            # Remove dictionary item
            result_vec.occup_basis.remove(index) # Note: we are not storing each basis state index
            prefactor = 1
            return (result_vec, prefactor)
        else:
            prefactor = np.sqrt(result_vec.occups[index])
            result_vec.occups[index] -= 1
            return (result_vec, prefactor)
    
    def creation(self,index):
        
        '''
        Method to carry out creation on a Fock ket index:
        Method doesn't modify original fock vector
        Returns:
                result_vec = Fock vector
                prefactor = creation factor
        '''
        assert (index >= 0 and index < self.M)
        
        result_vec = copy.deepcopy(self)
        prefactor = 1
        if (index not in self.occup_basis):   
            result_vec.occups[index] = 1
            result_vec.occup_basis.append(index)
        else:
            prefactor = np.sqrt(result_vec.occups[index] + 1)
            result_vec.occups[index] += 1
        
        return (result_vec, prefactor)
    
    
    def print_info(self):
        '''
        Print information on Fock vector
        '''
        print('Index: ', self.index)
        value = np.zeros(self.M)
        for i in self.occup_basis:
            value[i] = self.occups[i]
        print('State: ', value)
        return value
        
def test_init(N, M, occupations):
    '''
    Testing initialisation
    '''
    test_vector = fock_vector(N, M, occupations)
    print('Input')
    print('N %d M %d'%(N, M))
    print(occupations)
    print('---------------')
    print('Occupied bases: ', test_vector.occup_basis)
    print('Occupations:', test_vector.occups)
    
def test_annihilate(fock_vector, index):

    test_vector = fock_vector
    print('Input')
    test_vector.print_info()
    print('Annihilation on index ', index)
    print('Resulting Fock vector')
    print(test_vector.annihilate(index)[0].print_info())
    print('Prefactor: ', test_vector.annihilate(index)[1])
    print()
    print()
    
def test_create(fock_vector, index):
    test_vector = fock_vector
    print('Input')
    test_vector.print_info()
    print('Creation on index ', index)
    print('Resulting Fock vector')
    print(test_vector.creation(index)[0].print_info())
    print('Prefactor: ', test_vector.creation(index)[1])
    
# Simple Tests

#test_init(5, 3, [1,2,2])
#test_init(3,5, [1, 1, 0, 0, 1])

#test_annihilate(fock_vector(5, 3, [1,2,2]), 1)
#test_annihilate(fock_vector(3,5, [1, 1, 0, 0, 1]), 2)
#test_annihilate(fock_vector(3,5, [1, 0, 0, 0, 2]),2)
#test_create(fock_vector(3,2,[1,2]),2)
#test_create(fock_vector(3,2,[2,1]),2)
#test_create(fock_vector(3,2,[1,1]),2)
'''
vector = fock_vector(3, 2, [2, 1])
vector.print_info()
print('Annihilation on 0')
vector_1, prefactor_1 = vector.annihilate(0)
vector_1.print_info()
print(prefactor_1)
print('Creation on 0')
vector_2, prefactor_2 = vector_1.creation(0)
vector_2.print_info()
print(prefactor_2)
print('prefactor ', prefactor_1*prefactor_2)
'''
