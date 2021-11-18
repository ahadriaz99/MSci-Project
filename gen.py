# Functions to generate bosonic partitions (ie. repeated combinations
# Written by Lajos Palanki, used with permission from author

import numpy as np
import math

def get_sum(n, m):
    if (n == 0):
        return 1
    if (m == 0):
        return 0
    if (m == 1):
        return n
    
    s = 0
    for i in range(n):
        s += get_sum(n-i,m-1) + 1
    return s

def get_partial_sum(n, m, mx):
    s = 0
    if (n <= 0):
        return 1
    if (m <= 0):
        return 0
    
    for i in range(mx):
        s += get_sum(n-i,m-1) + 1
    return s

def get_highest_occup(n, id):
    m = 0
    while(id > get_sum(n, m)):
        m += 1
    
    return m

def get_highest_occup_number(n, id, m = None):
    if(id == 0 and (m == None or m == 0)):
        return n
    if(id == 0 and m != 0):
        return 0
    max_occup = get_highest_occup(n,id)
    offset_id = id - get_sum(n,max_occup-1) - 1
    occup = 0
    if(m == None or max_occup == m):
        occup = 1
    if(m == None):
        m = max_occup
    if(n != 1):
        return occup + get_highest_occup_number(n-1, offset_id, m)
    return occup

def get_basis(n,id):
    mn = get_highest_occup_number(n,id)
    mid = get_highest_occup(n,id)
    ret = np.zeros(mid + 1)
    ret[mid] = mn
    if (n == 1):
        return ret
    if (mid == 0):
        return ret
    mxid = 0

    while id >= get_partial_sum(n,mid,mxid+1):
        mxid += 1

    nid = id - get_partial_sum(n,mid,mxid)
    if n - mn != 0:
        r = get_basis(n - mn, nid)
        for i,val in enumerate(r):
            ret[i] = val

    return ret



if __name__ == "__main__":

    n = 3
    m = 3
    counter = 0
    for id in range(get_sum(n,m)+1):
        print(get_basis(n,id))
        counter += 1
    m += 1
    print(counter, int(math.factorial(n+m-1)/math.factorial(m-1)/math.factorial(n)))
