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

def ruleAscLen(n, l):
# Reproduced from [https://jeromekelleher.net/category/combinatorics.html] and
# Reproduced from
#[https://math.stackexchange.com/questions/18659/algorithm-for-generating-integer-partitions]
    a = [0 for i in range(n + 1)]
    k = 1
    a[0] = 0
    a[1] = n
    while k != 0:
        x = a[k - 1] + 1
        y = a[k] - 1
        k -= 1
        while x <= y and k < l - 1:
            a[k] = x
            y -= x
            k += 1
        a[k] = x + y
        yield a[:k + 1]

        
N = 10
M = 5
L_max = N*(M-1) # Check entirety of Hilbert space

# Basis generation
configs.disc_config_very_fast(N, M, L_max)

fock_size_H = 0
fock_size_H1 = 0
fock_size_partition = 0 # Generating fock space from partitions
for L in range(L_max+1):
    
    print('Simulation parameters')
    print('Disc geometry')
    print('N ', N, ' M ', M, ' L ', L)
    fast_t0 = timeit.default_timer()

    H = disc_Hamiltonian_fast(N=N,M=M,L=L)
    H.generate_basis()
    # Add fock size of subspace
    fock_size_H += H.fock_size
    #H.construct_Hamiltonian_fast()

    fast_t1 = timeit.default_timer()

    reg_t0 = timeit.default_timer()

    H1 = disc_Hamiltonian(N=N,M=M,L=L)
    H1.generate_basis()
    #H1.construct_Hamiltonian()
    # Add fock size of subspace

    fock_size_H1 += H1.fock_size

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
    
    # Count the number of limited partitions
    if (L > 1):
        count = 0
        for part in ruleAscLen(L-1, N):
            fock_size_partition += 1
            count += 1
        
        assert H.fock_size == H1.fock_size
        if (count != H.fock_size):
            print('Generated partition dimension is not equivalent to Fock size.')
        print('Reference Fock size ', H.fock_size)
        print('Partition dimension ', count)
        print('Basis')
        H.show_basis()
        print('Partitions')
        for part in ruleAscLen(L-1, N):
            print(part)


if (fock_size_H != math.factorial(N+M-1)/math.factorial(N)/math.factorial(M-1) or \
    fock_size_H1 != math.factorial(N+M-1)/math.factorial(N)/math.factorial(M-1)):
    print('Basis dimensions are wrong!')
    print('Expected Fock space size: ', int(math.factorial(N+M-1)/math.factorial(N)/math.factorial(M-1)))
    print('H size ', fock_size_H)
    print('H1 size ', fock_size_H1)
    assert 0

print('Basis set size matches analytical Fock size.')

# partition cannot generate Bose condensed basis state, add one to fock size
if (fock_size_partition+1 != math.factorial(N+M-1)/math.factorial(N)/math.factorial(M-1)):
   
    print('Partition dimensions are wrong!')
    print('Expected Fock space size: ', int(math.factorial(N+M-1)/math.factorial(N)/math.factorial(M-1)))
    print('Partition size ', fock_size_partition)
    
