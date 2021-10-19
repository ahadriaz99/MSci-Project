# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 15:20:35 2021

@author: Owner
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

#Test case
#configurations(2,2)
#configurations(3,2)