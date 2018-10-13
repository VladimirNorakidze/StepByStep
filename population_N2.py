# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 17:53:38 2018

@author: Владимир
"""
import constants as cnst
import numpy as np

def N_rot_trans(B, diss):
    return int(0.5 * (-1 + np.sqrt(1 + 4 * float(diss)/cnst.B)))

N = N_rot_trans(cnst.B, cnst.diss_energy_N2)+50
g = np.ndarray(N+5)
x_J = np.ndarray(N+5)
eps_J = np.ndarray(N+5)
eps_th_ex = np.ndarray(N+5)
g[0] = 6
g[1] = 9
x_J[0] = (2./9)*(cnst.B*11600./cnst.T)*6
x_J[1] = (2./9)*(cnst.B*11600./cnst.T)*9*np.exp(-eps_J[1])
for J in range(N):
    eps_J[J] = cnst.B*J*(J+1)*11600./cnst.T #Безразмерная энергия уровня J
for J in range(N):
    eps_th_ex[J] = eps_J[J+2]-eps_J[J]
for J in range(2, N):
    if (J%2.==0):
        g[J]   = (2*J+1)*6 #Взято из статьи про HITRAN
        x_J[J] = (2./9)*(cnst.B*11600./cnst.T)*g[J]*np.exp(-eps_J[J]) #заселенность J уровня
    else:
        g[J]   = (2*J+1)*3
        x_J[J] = (2./9)*(cnst.B*11600./cnst.T)*g[J]*np.exp(-eps_J[J])      
