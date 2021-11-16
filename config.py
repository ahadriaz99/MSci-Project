# -*- coding: utf-8 -*-
"""
Function to generated repeated combinatrics of bosonic configuration

19/10/21

by Marcell Dorian Kovacs and Ahad Riaz
"""
import numpy as np
from itertools import combinations
from Fock_vector import fock_vector
from os.path import exists
from os import mkdir
import os
from os import system

def configurations(N, M):
    ''' the function is used to implement a algorithm to find all the possible 
    configurations for a given number of particles (N) and number of states (M)
    '''
    valid_configs = list()
    for comb in combinations(range(N + M - 1), M - 1):
        valid_configs.append(np.flip([b - a - 1 for a, b in zip((-1,) + comb, comb + (N + M - 1,))]))
            
    valid_configs = list(map(tuple,valid_configs))
    return valid_configs

def disc_config(N,M,L):
    
    #print('configurations', config_input)
    for i in range(0,L+1):
        if exists('Disc_Configurations_N%dM%dL%d.txt'%(N, M, i)):
            #print('Configurations have already been generated for N %d M %d L %d '%(N, M, i))
            if (i == L):
                return
            continue
        with open('Disc_Configurations_N%dM%dL%d.txt'%(N, M, i), 'a') as f:
            f.write('N,   M,    L,    Basis \n')
    f.close()
    config_input = np.array(configurations(N, M)) # Calculate repeated combinations
    
    for j in range(len(config_input)):
        assert len(config_input[j]) == M
        vector = fock_vector(N, M, config_input[j])
        ang_mom = vector.ang_mom()
        if ang_mom > L:
            continue
        else:
            with open('Disc_Configurations_N%dM%dL%d.txt'%(N, M, ang_mom), 'a') as f:
                f.write('%d   %d   %d     '%(N,M,ang_mom))
                for num in config_input[j]:
                    f.write(str(num)+' ')
                f.write('\n')
    f.close()
    
#Test case
#disc_config(4, 7, 12)
#configurations(2,2)
#configurations(3,2)
