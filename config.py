# -*- coding: utf-8 -*-
"""
Function to generated repeated combinatrics of bosonic configuration

19/10/21

by Marcell Dorian Kovacs and Ahad Riaz
"""
import numpy as np
from itertools import combinations

def configurations(N, M):
    ''' the function is used to implement a algorithm to find all the possible 
    configurations for a given number of particles (N) and number of states (M)
    '''
    valid_configs = list()
    for comb in combinations(range(N + M - 1), M - 1):
        valid_configs.append(np.flip([b - a - 1 for a, b in zip((-1,) + comb, comb + (N + M - 1,))]))
            
    valid_configs = list(map(tuple,valid_configs))
    return valid_configs

def config_ang(N,M,L):
    
    config_input = np.array(configurations(N, M)) # Calculate repeated combinations
    #print('configurations', config_input)
        
    index = 0
    basis = []
    for i in range(len(config_input)):
        #print(i)
        assert len(config_input[i]) == M # Check correct input format
        vector = fock_vector(N, M, config_input[i])
        # Only add restricted ang. mom. bases to the Fock spaces
        if(vector.ang_mom() == L):
            basis.append(fock_vector(N, M, config_input[i],index=index)) # Create fock vectors
            index += 1
          
    
    for basis in basis:
        return(basis.print_info())
    
    
    with open('Restricted_Bases_N%dM%d.txt'%(N, M), 'a') as f:
            f.write('&N= %d, M= %d, L=%d, \n&END\n'%(self.M, self.N))
#Test case
#configurations(2,2)
#configurations(3,2)