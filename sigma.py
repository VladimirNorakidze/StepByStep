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

def r_excitation(energy, J0, sw=None):
    sigma = 0.
    treshold = threshold_energy['rot'][J0]
    if (energy > treshold): #Если меньше, то сечение равно нулю
        sigma = (8./15)*np.pi*(const.Q*const.a0)**2*np.sqrt(1-treshold/energy)*( ( (J0+1)*(J0+2) )/( (2*J0+1)*(2*J0+3) ) )
    return sigma/const.sigma0

def r_deexcitation(energy, J0, sw=None):
    sigma = 0.
    treshold = threshold_energy['rot'][J0]
    if energy > 0.:
        sigma = (8./15)*np.pi*(const.Q*const.a0)**2*np.sqrt(1+treshold/energy)*( ( (J0+1)*(J0+2) )/( (2*J0+5)*(2*J0+3) ) )
    return sigma/const.sigma0
