# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 16:59:33 2018

@author: Владимир
"""
import sigma
import json

"""Загрузка пороговых энергий неупругих переходов молекулы азота"""
with open("threshold_energy.txt","r") as file:
    threshold_energy = json.load(file)
    
def rotational_term(eps, J0, x_J, f0, f0D, func):
    treshold = threshold_energy['rot'][J0]
    if (f0D[10] == 1):
        F_ex   = -x_J[J0]   * f0 * (eps*sigma.r_excitation(eps,J0)   - \
                     (eps+treshold)*sigma.r_excitation(eps+treshold,J0))
        F_deex = -x_J[J0+2] * f0 * (eps*sigma.r_deexcitation(eps,J0) - \
                     (eps-treshold)*sigma.r_deexcitation(eps-treshold,J0))
    else:
        f1 = func(eps+treshold)
        if eps-treshold <= 0:
            f2 = 0
        else:
            f2 = func(eps-treshold)
        F_ex   = -x_J[J0]   * (eps*sigma.r_excitation(eps,J0)*func(eps)   - \
                     (eps+treshold)*sigma.r_excitation(eps+treshold,J0)*f1)
        F_deex = -x_J[J0+2] * (eps*sigma.r_deexcitation(eps,J0)*func(eps) - \
                     (eps-treshold)*sigma.r_deexcitation(eps-treshold,J0)*f2)
    return [F_ex, F_deex] 

def vibrational_term(eps, f0, v, k, f0D, func):
    threshold = threshold_energy['vibr'][v][k]
    if (f0D[0] == 1):
        F_ex   = -0.1   * f0 * (eps*sigma.v_excitation_vk(eps, v, k)   - (eps+threshold)*sigma.v_excitation_vk(eps+threshold, v, k))
        F_deex = -0.1 * f0 * (eps*sigma.v_deexcitation_vk(eps, v, k) - (eps-threshold)*sigma.v_deexcitation_vk(eps-threshold,v,k))
    else:
        f1 = func(eps + threshold)
        if eps-threshold <= 0:
            f2 = 0
        else:
            f2 = func(eps - threshold)
        F_ex   = -0.1   * (eps*sigma.v_excitation_vk(eps,v,k)*func(eps)   - (eps+threshold)*sigma.v_excitation_vk(eps+threshold,v,k)*f1)
        F_deex = -0.1 * (eps*sigma.v_deexcitation_vk(eps,v,k)*func(eps) - (eps-threshold)*sigma.v_deexcitation_vk(eps-threshold,v,k)*f2)
    return [F_ex, F_deex] 