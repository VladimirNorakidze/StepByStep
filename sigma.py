# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 16:27:57 2018

@author: Владимир
"""
import numpy as np
from scipy import interpolate
import json

import constants as const

"""Загрузка пороговых энергий неупругих переходов молекулы азота"""
with open("threshold_energy.txt","r",encoding="utf-8") as file:
    threshold_energy = json.load(file)

"""Значения сечения упругого взаимодействия и энергии из Итикавы:"""
eps  = np.array([0,0.55,0.7,0.9,1.0,1.5,2.0,2.2,2.35,2.5,2.7,
                 3.0,4.0,5.0,6.0,8.0,10,15,20,25,30,40,50,60,80,
                 100,120,150,200,250,300,400,500,600,800,1000,1e17])*(11600./const.T)
    
sig  = np.array([0,8.39,9.03,9.62,9.83,10.53,17.93,19.5,20.5,21.0,17.5,
                 15.0,11.6,10.75,10.6,10.6,11.4,11.8,11.15,10.25,9.65,8.85,
                 8.2,7.4,6.25,5.6,4.9,4.2,3.5,3.0,2.65,2.15,1.85,1.6,1.25,1.0,0])/(const.sigma0*1.0e20)
    
"""Функция сечения упругого взаимодействия:"""
elastic = interpolate.UnivariateSpline(eps, sig, s=0)

def r_excitation(energy, J0):
    sigma = 0.
    treshold = threshold_energy['rot'][J0]
    if (energy > treshold): #Если меньше, то сечение равно нулю
        sigma = (8./15)*np.pi*(const.Q*const.a0)**2*np.sqrt(1-treshold/energy)*( ( (J0+1)*(J0+2) )/( (2*J0+1)*(2*J0+3) ) )
    return sigma/const.sigma0

def r_deexcitation(energy, J0):
    sigma = 0.
    treshold = threshold_energy['rot'][J0]
    if energy > 0.:
        sigma = (8./15)*np.pi*(const.Q*const.a0)**2*np.sqrt(1+treshold/energy)*( ( (J0+1)*(J0+2) )/( (2*J0+5)*(2*J0+3) ) )
    return sigma/const.sigma0

"""Колебательное возбуждение"""

a_coefs = np.array([1.83e10, 10., 9.41e-4, 0.42, 1.24e9, 10., 2.087e-3, 7.98, 1.37e-2, 9.2, 1.94e-2, 6.9])

def f1 (t, c1, c2):
    return 1e-20 * c1 * (t/const.E_R)**c2

def f2 (t, c1, c2, c3, c4):
    return f1(t, c1, c2) / ( 1 + (t/c3)**(c2+c4) )

def v_excitation_01(energy): #Функция сечения возбуждения колебательных степеней для перехода 0->1:
    s1 = f2((energy - threshold_energy["vibr"][0][1])*(873./11600)*1e-3, a_coefs[0], a_coefs[1], a_coefs[2], a_coefs[3])
    s2 = f2((energy - threshold_energy["vibr"][0][1])*(873./11600)*1e-3, a_coefs[4], a_coefs[5], a_coefs[6], a_coefs[7])
    s3 = f2((energy - threshold_energy["vibr"][0][1])*(873./11600)*1e-3, a_coefs[8], a_coefs[9], a_coefs[10], a_coefs[11])
    return (s1 + s2 + s3) / const.sigma0

factors_vcs = np.array([0, 1, 0.49, 0.45, 0.2, 0.13, 0.54e-1, 0.19e-1, 0.73e-2,\
                        0.21e-2, 0.5e-3, 0.98e-4, 0.17e-4, 0.28e-5, 0.67e-6, \
                        0.14e-6, 0.25e-7, 0.39e-8, 0.56e-9, 0.86e-10])

def v_excitation_0v(energy, v):
    if energy > threshold_energy['vibr'][0][v]:
        return (factors_vcs[v]) * v_excitation_01(energy)
    else:
        return 0
    
def v_excitation_vk(energy, v, k):
    sigma = 0
    if energy > threshold_energy['vibr'][v][k]:
        eps = energy + threshold_energy['vibr'][0][v]
        sigma_0v = v_excitation_0v(eps, v)
        sigma_0k = v_excitation_0v(eps, k)
        if v == 0:
            sigma_0v = elastic(eps)
        elif k == 0:
            sigma_0k = elastic(eps)
        sigma = (eps/energy) * sigma_0v * sigma_0k / elastic(eps)
    return sigma
#sigma_vibr_ex = interpolate.UnivariateSpline(eps, sig, s=0)

def v_deexcitation_vk(energy,v,k): #Функция для перехода k->v
    eps = energy + threshold_energy['vibr'][v][k]
    sigma = v_excitation_vk(eps, v, k) * eps/energy
    return sigma