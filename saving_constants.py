# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 14:35:31 2018

@author: Владимир
"""

def save_constants(U):
    with open("constants.py", 'w') as f:
        f.write("""import numpy as np
e = 1.6*10**(-19)
m = 9.1*10**(-31)
M0 = 1.67*10**(-27)
k = 1.38e-23
c = 3.0e8
Q = -1.13
alpha = 3.13
Planck_const = 1.0545718e-34
a0 = 137.*Planck_const/(m*c)
B = 2.49e-4
diss_energy_N2 = 9.73
A = 28
U = {}
r = 0.9
E = U/r
M = A*M0
delta = 2.*m/M
p = 1.5*133
T = 873.
n0 = 1.0e16
sigma0 = np.pi*a0**2
Part_Kn = np.sqrt(1./3)*e*E/( (k*T) * sigma0)
E_R = 13.61
N_all  = p/(k*T)""".format(U))