# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 14:09:25 2021

@author: Owner
"""

from Hamiltonian import Hamiltonian

#H = ham.Hamiltonian(2,3)


#H.generate_basis()
#H.show_basis()
#H.construct_Hamiltonian()
#H.print_matrix(H.many_body_H)

#H.show_basis()
#H.print_matrix(H.basis_overlap())

#evalues, evecs = H.diagonalise()
#print('Hamiltonian eigenvalues [V0]')
#print(evalues)
#print('Ground state energy [V0] ', H.e_ground)

#H.check_sign_problem()

def Dirac_Delta(a, b):
    '''
    Simple delta function
    '''
    if (a == b):
        return 1
    else:
        return 0
    
class sub_Hamiltonian(Hamiltonian):
    
    def __init__(self,N,M):
        super().__init__(N,M)
          
    def matrix_overlap_disc(self, i, j, k, l):
        '''
        Construct many-body matrix elements for disc Hamiltonian
        '''
        V0 = 1
        return Dirac_Delta(i+j, k+l)*V0*math.factorial(i+j)/2**(i+j)/\
               np.sqrt(math.factorial(i)*math.factorial(j)*math.factorial(k)*math.factorial(l))
        
    def H_element(self, basis1, basis2):
        '''
        Calculate matrix element between 2 many-body basis states
        '''
        assert basis1 in self.basis and basis2 in self.basis
        #print('Basis 1')
        #basis1.print_info()
        #print('Basis 2')
        #basis2.print_info()
        element = 0 # Matrix element
        # Loop over all possible single-particle state indices
        # NEEDS OPTIMISATION
        for i in range(self.M):
            for j in range(self.M):
                for k in range(self.M):
                    for l in range(self.M):
                    
                        #print('Operator indices', i, j, k, l)
                        #print('BRA BASIS')
                        #basis1.print_info()
                        #print('KET BASIS')
                        #basis2.print_info()
                        
                        # Calculate scalar integral overlaps
                        matrix_overlap = self.matrix_overlap_disc(i, j, k, l)
                        if (matrix_overlap != 0):
                            # Apply annihlation on basis 1 BRA -> (i, j) must be reversed
                            new_basis1, total_prefactor_1 = self.annihilate(basis1, j, i)
                            # Apply annihilation on basis 2 KET
                            new_basis2, total_prefactor_2 = self.annihilate(basis2, k, l)
                            #print('TRANSFORMED BASIS1')
                            #new_basis1.print_info()
                            #print('TRANSFORMED BASIS2')
                            #new_basis2.print_info()
                            #print('Prefactor ', total_prefactor_1*total_prefactor_2)
                            # Assemble matrix element contribution
                            element += 0.5*matrix_overlap*self.overlap(new_basis1, new_basis2)*total_prefactor_1*total_prefactor_2
                            #print('Overlap: ', self.overlap(new_basis1, new_basis2))
        return element

H = sub_Hamiltonian(4,2)


H.generate_basis()
H.show_basis()
H.construct_Hamiltonian()
H.print_matrix(H.many_body_H)

#H.show_basis()
#H.print_matrix(H.basis_overlap())