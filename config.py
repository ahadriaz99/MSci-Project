# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 15:20:35 2021

@author: Owner
"""

def configurations(N, M):
    ''' the function is used to implement a algorithm to find all the possible 
    configurations for a given number of particles (N) and number of states (M)
    '''
    configsFound = list() #initiate the valid configuration list

    def findValidConfigs(particlesNeeded, totalStates, currentState,
                         currentConfig):
        '''

        Parameters
        ----------
        particlesNeeded : the number of particles
        totalStates : the total number of states
        currentState : the current level of the state
        currentConfig : the most recent configuration


        '''
        if currentState > totalStates:
            if particlesNeeded == 0:
                configsFound.append(currentConfig)
                return
            else:
                return

        for i in range(particlesNeeded + 1):
            findValidConfigs(particlesNeeded - i, totalStates,
                             currentState + 1, 
                             currentConfig + (i,))

    findValidConfigs(N, M, 1, tuple())
    return configsFound

#Test case
#configurations(2,2)
#configurations(3,2)