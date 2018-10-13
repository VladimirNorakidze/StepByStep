# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 16:59:33 2018

@author: Владимир
"""
import sigma

def rotational_term(eps, J0, x_J, f0, f0D, treshold, func):
    if (f0D[10] == 1):
        F_ex   = -x_J[J0]   * f0 * (eps*sigma.r_excitation(eps,J0,treshold)   - \
                     (eps+treshold)*sigma.r_excitation(eps+treshold,J0,treshold))
        F_deex = -x_J[J0+2] * f0 * (eps*sigma.r_deexcitation(eps,J0,treshold) - \
                     (eps-treshold)*sigma.r_deexcitation(eps-treshold,J0,treshold))
    else:
        f1 = func(eps+treshold)
        if eps-treshold <= 0:
            f2 = 0
        else:
            f2 = func(eps-treshold)
        F_ex   = -x_J[J0]   * (eps*sigma.r_excitation(eps,J0,treshold)*func(eps)   - \
                     (eps+treshold)*sigma.r_excitation(eps+treshold,J0,treshold)*f1)
        F_deex = -x_J[J0+2] * (eps*sigma.r_deexcitation(eps,J0,treshold)*func(eps) - \
                     (eps-treshold)*sigma.r_deexcitation(eps-treshold,J0,treshold)*f2)
    return [F_ex, F_deex] 