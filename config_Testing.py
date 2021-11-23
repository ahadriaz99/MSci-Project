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
import timeit

from Fock_vector import fock_vector
from Disc_Hamiltonian import disc_Hamiltonian
from Disc import disc_Hamiltonian_fast
import Ryser_Algorithm as ryser
import config as configs
from numpy import linalg as la
    
def check_matrix(tolerance, M1, M2):
    assert M1.shape == M2.shape
    for i in range(len(M1)):
        for j in range(i, len(M1)):
            if (abs(M1[i, j] - M1[j, i]) >\
                   tolerance or \
                abs(M2[i, j] - M2[j, i]) >\
                   tolerance or \
                abs(M1[i, j] - M2[i, j]) >\
                   tolerance):
                print('Element mismatch!')
                print(i, j, M1[i, j], M1[j, i])
                print(i, j, M2[i, j], M2[j, i])
                assert 0


        
N = 5
M = N*(N-1)
L_max = M

# Basis generation
configs.disc_config_very_fast(N, M, L_max)

for L in range(L_max):
    
    print('Simulation parameters')
    print('Disc geometry')
    print('N ', N, ' M ', M, ' L ', L)
    fast_t0 = timeit.default_timer()

    H = disc_Hamiltonian_fast(N=N,M=M,L=L)
    H.generate_basis()
    #H.construct_Hamiltonian_fast()

    fast_t1 = timeit.default_timer()

    reg_t0 = timeit.default_timer()

    H1 = disc_Hamiltonian(N=N,M=M,L=L)
    H1.generate_basis()
    #H1.construct_Hamiltonian()

    reg_t1 = timeit.default_timer()
    



    print('Time for fast routine [s]: ', fast_t1 - fast_t0)
    print('Time for regular routine [s]: ', reg_t1 - reg_t0)
    print('Speedup ', ((reg_t1 - reg_t0)/(fast_t1 - fast_t0)))

    print('Running basis comparison...')
    
    equiv = np.zeros(len(H.basis))
    for i in range(len(H.basis)):
        for basis2 in H1.basis:
            if (H.basis[i].occups == basis2.occups):
                equiv[i] = 1
                continue
        if equiv[i] == 0:
            print('Fast routine basis not contained.')
            H.basis[i].print_info()
            print('\n\nFast routine basis\n\n')
            H.show_basis()
            print('\n\nOriginal routine basis\n\n')
            H1.show_basis()
            assert 0
            
    print('Basis sets are equivalent up to reordering.')
