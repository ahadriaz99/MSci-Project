# -*- coding: utf-8 -*-
"""
Function to generated repeated combinatrics of bosonic configuration

19/10/21

by Marcell Dorian Kovacs and Ahad Riaz
"""
import numpy as np
import gen as generator
import math
from itertools import combinations
from Fock_vector import fock_vector
from os.path import exists
from os import mkdir
import os
from os import system

from const_l_gen import gen_bases

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
    
def disc_config_very_fast(N, M, L):
    print('N, M, L ', N, M, L )
    m = M-1
    counter = 0

    for i in range(0,L+1):
        if exists('Disc_Configurations_N%dM%dL%d.txt'%(N, M, i)):
            #print('Configurations have already been generated for N %d M %d L %d '%(N, M, i))
            if (i == L):
                return
            continue
        with open('Disc_Configurations_N%dM%dL%d.txt'%(N, M, i), 'a') as f:
            f.write('N,   M,    L,    Basis \n')
        f.close()
        
        bases = gen_bases(N, i)

        with open('Disc_Configurations_N%dM%dL%d.txt'%(N, M, i), 'a') as f:
            for basis in bases:
                print('Raw generated basis: ',basis)
                if len(basis) > M:
                    #basis = basis[:M]
                   # assert len(basis) == M

                    if np.array(basis[M:]).any() != 0:
                        continue
                    else:
                        basis = basis[:M]
                    
                f.write(str(N)+' '+str(M)+' '+str(i)+'    ')
                if len(basis) < M:
                    basis = np.array(np.concatenate((np.array(basis),np.zeros(abs(len(basis)-M)))))
                    basis = basis.astype(int)
                    assert len(basis) == M

                for j in range(len(basis)):
                    f.write(str(basis[j])+' ')
                f.write('\n')
                
                #print('Formatted basis: ',basis)
                counter += 1
        f.close()
                
        print('Basis generation progress [%] ', (i/L)*100)
            
            
    
def disc_config_fast(N, M, L):
    n = N
    m = M-1
    counter = 0
    
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
    
    subspace_size = np.zeros(L+1)
    for id in range(generator.get_sum(n,m)+1):
        basis = np.array(generator.get_basis(n,id))
        if len(basis) < M:
            basis = np.array(np.concatenate((np.array(basis), np.zeros(abs(len(basis)-M)))))
        basis = basis.astype(int)

        vector = fock_vector(N, M, basis)
        ang_mom = vector.ang_mom()
        if ang_mom > L:
            if (counter % 1000 == 0):
                size = int(math.factorial(N+M-1)/math.factorial(M-1)/math.factorial(N))
                print('Basis generation progress [%] ', (counter/size)*100)
            counter += 1
            continue
        else:
            with open('Disc_Configurations_N%dM%dL%d.txt'%(N, M, ang_mom), 'a') as f:
                f.write('%d   %d   %d     '%(N,M,ang_mom))
                for num in basis:
                    f.write(str(num)+' ')
                f.write('\n')
            subspace_size[ang_mom] += 1
        
        if (counter % 10000 == 0):
            size = int(math.factorial(N+M-1)/math.factorial(M-1)/math.factorial(N))
            print('Basis generation progress [%] ', (counter/size)*100)
        counter += 1
    
    with open('Stats_Disc_Configurations_N%dM%dL%d.txt'%(N, M, L), 'a') as f:
                f.write('N = %d M =  %d L_max =  %d     \n'%(N,M,L))
                for index in range(len(subspace_size)):
                    f.write('L= '+str(index)+' subspace size: '+str(subspace_size[index])+' \n')
                f.write('\n')
                
        
    assert counter == int(math.factorial(N+M-1)/math.factorial(M-1)/math.factorial(N))

def sphere_config_fast(N, M, L, S):
    assert M == 2*S + 1
    S0 = S
    n = N
    m = M-1
    counter = 0

    #print('configurations', config_input)
    for i in range(0,L+1):
        if exists('Sphere_Configurations_N%dM%dL%dS%d.txt'%(N, M, i, S0)):
            #print('Configurations have already been generated for N %d M %d L %d '%(N, M, i))
            if (i == L):
                return
            continue
        with open('Sphere_Configurations_N%dM%dL%dS%d.txt'%(N, M, i, S0), 'a') as f:
            f.write('N,   M,    L,    S,    Basis \n')
    f.close()
    
    subspace_size = np.zeros(L+1)
    for id in range(generator.get_sum(n,m)+1):
        basis = np.array(generator.get_basis(n,id))
        print(basis)
        if len(basis) != M:
            basis = np.array(np.concatenate((np.array(basis), np.zeros(abs(len(basis)-M)))))
        basis = basis.astype(int)
        print(basis)
        vector = fock_vector(N, M, basis, S = S0)
        ang_mom = vector.ang_mom()
        if ang_mom > L or ang_mom < 0:
            if (counter % 1000 == 0):
                size = int(math.factorial(N+M-1)/math.factorial(M-1)/math.factorial(N))
                print('Basis generation progress [%] ', (counter/size)*100)
            counter += 1
            continue
        
        else:
            with open('Sphere_Configurations_N%dM%dL%dS%d.txt'%(N, M, ang_mom, S0), 'a') as f:
                f.write('%d   %d   %d   %d  '%(N,M,ang_mom, S0))
                for num in basis:
                    f.write(str(num)+' ')
                f.write('\n')
            subspace_size[ang_mom] += 1
        
        if (counter % 10000 == 0):
            size = int(math.factorial(N+M-1)/math.factorial(M-1)/math.factorial(N))
            print('Basis generation progress [%] ', (counter/size)*100)
        counter += 1
    
    with open('Stats_Sphere_Configurations_N%dM%dL%dS%d.txt'%(N, M, L, S0), 'a') as f:
                f.write('N = %d M =  %d L_max =  %d  S = %d   \n'%(N,M,L,S0))
                for index in range(len(subspace_size)):
                    f.write('L= '+str(index)+' subspace size: '+str(subspace_size[index])+' \n')
                f.write('\n')
        
    assert counter == int(math.factorial(N+M-1)/math.factorial(M-1)/math.factorial(N))

def sphere_config_very_fast(N, M, L, S):
    M = 2*S + 1 
    S0 = S
    m = M-1
    counter = 0
    count = 0
    for i in range(0,L+1):
        if exists('Sphere_Configurations_N%dM%dL%dS%d.txt'%(N, M, i, S0)):
            #print('Configurations have already been generated for N %d M %d L %d '%(N, M, i))
            if (i == L):
                return
            continue
        with open('Sphere_Configurations_N%dM%dL%dS%d.txt'%(N, M, i, S0), 'a') as f:
            f.write('N,   M,    L,    S,    Basis \n')
        f.close()
        
        bases = gen_bases(N, i)
        
        for basis in bases:
            print (basis)
            if len(basis) > M:
                     #basis = basis[:M]
                    # assert len(basis) == M
 
                if np.array(basis[M:]).any() != 0:
                    continue
                else:
                    basis = basis[:M]
            print(basis)
            if len(basis) < M:
                basis = np.array(np.concatenate((np.array(basis), np.zeros(abs(len(basis)-M)))))
                basis = basis.astype(int)
                assert len(basis) == M
            print(basis)
            
            vector = fock_vector(N,M,basis, S = S0)
            ang_mom = vector.ang_mom()
            print(ang_mom)
            count += 1
            
            if ang_mom > L or ang_mom < 0:
                if (counter % 1000 == 0):
                    size = int(math.factorial(N+M-1)/math.factorial(M-1)/math.factorial(N))
                    print('Basis generation progress [%] ', (counter/size)*100)
                counter += 1
                continue
       
            else:
               with open('Sphere_Configurations_N%dM%dL%dS%d.txt'%(N, M, ang_mom, S0), 'a') as f:
                   f.write('%d   %d   %d   %d  '%(N,M,ang_mom, S0))
                   for num in basis:
                       f.write(str(num)+' ')
                   f.write('\n')
            #print(basis)
            ##vector = fock_vector(N,M,basis, S =S0)
            #ang_mom = vector.ang_mom()
            #print(basis, ang_mom)
            
                
#sphere_config_fast(10,5,10,5)
#disc_config_very_fast(3,5 , 5)
#Test case
#sphere_config_fast(3, 5, 5,2)
sphere_config_very_fast(3, 5, 15,2)
#configurations(2,2)
#configurations(3,2)
