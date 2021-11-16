# -*- coding: utf-8 -*-
"""
Function to generated repeated combinatrics of bosonic configuration

19/10/21

by Marcell Dorian Kovacs and Ahad Riaz
"""
import numpy as np
from itertools import combinations
from Fock_vector import fock_vector

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
    
    for i in range(0,L+1):
        with open('Configurations_N%dM%dL%d.txt'%(N, M, i), 'a') as f:
            f.write('N,   M,    L,    Basis \n')
        f.close()
    
    for j in range(len(config_input)):
        assert len(config_input[i]) == M
        vector = fock_vector(N, M, config_input[i])
        ang_mom = vector.ang_mom()
        value = np.zeros(M)
        for k in vector.occup_basis:
            value[k] = vector.occups[k]
        with open('Configurations_N%dM%dL%d.txt'%(N, M, ang_mom), 'a') as f:
             f.write('%d,     %d,    %d, %s \n'%(N,M,ang_mom,value))
        
#Test case
config_ang(2,3,3)
#configurations(2,2)
#configurations(3,2)