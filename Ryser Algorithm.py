# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 01:25:29 2021

@author: Owner
"""

from random import *
from time import time
from itertools import *

def perm(a): # for comparison
    if len(a) == 1:
        return a[0][0]
    elif len(a) == 2:
        return a[0][0]*a[1][1]+a[1][0]*a[0][1]
    else:
        tsum = 0
        for i in range(len(a)):
            transposed = [zip(*a)[j] for j in range(len(a)) if j != i]
            tsum += a[0][i] * perm(zip(*transposed)[1:])
        return tsum

def perm_ryser(a): # Ryser's formula, using matrix entries
    maxn = len(a)
    n_list = range(1,maxn+1)
    s_list = chain.from_iterable(combinations(n_list,i) for i in range(maxn+1))
    total = 0
    for st in s_list:
        stotal = (-1)**len(st)
        for i in range(maxn):
            stotal *= sum(a[i][j-1] for j in st)
        total += stotal
    return total*((-1)**maxn)


def genmatrix(d):
    mat = []
    for x in range(d):
        row = []
        for y in range(d):
            row.append([-5,5][randrange(0,2)])
        mat.append(row)
    return mat

