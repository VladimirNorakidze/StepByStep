# -*- coding: utf-8 -*-
"""
Редактор Spyder
Вид уравнения:сложный
Имеющиеся проблемы: 
    - при малых значениях U (напряжение) получается уравнение колебаний, чего быть не должно. 
        
"""
import numpy as np
from scipy.integrate import odeint
from scipy import interpolate
from scipy.integrate import simps
import matplotlib.pyplot as plt
from scipy.misc import derivative
import copy
import json
import time
import saving_constants as sc

start_time = time.time()
U = 1. # [В]
sc.save_constants(U)

import constants as c
import sigma
import inelastic_term as it
import population_N2 as pN2

"""Загрузка пороговых энергий неупругих переходов молекулы азота"""
with open("threshold_energy.txt","r",encoding="utf-8") as file:
    threshold_energy = json.load(file)
    
N  = 13000 #число точек
N_for_dif_case = 1.3e1
electron_energy = np.linspace( 1.0e-4, N_for_dif_case, N) # значения энергии
switch_vibr = 0
dt   = electron_energy[10] - electron_energy[9] #Шаг по энергии 
beta = 1.0e8
"""Расчет числа молекул. Они подчиняются закону распределения Максвелла"""
N_all  = c.p/(c.k*c.T) #м^-3

"""Функция для решения дифференцияального уравнения относительно f0:"""    
def search_f0(y,t):
    res   = np.ndarray((2))
    f0, p = y
    total_sigma = 0
    derivate_term = t*(2*sigma.elastic(t) + t*derivative(sigma.elastic,t))
    der_total_sigma_term = derivative(sigma.elastic,t)
    F_ex   = 0
    F_deex = 0
    V_Ex   = 0
    V_Deex = 0
    for J0 in range(pN2.N_rot_trans(c.B, c.diss_energy_N2)):
        Ex, Deex = it.rotational_term(t, J0, pN2.x_J, f0, f0D, func)
        F_ex   += Ex
        F_deex += Deex
        total_sigma          += pN2.x_J[J0] * sigma.r_excitation(t,J0) + pN2.x_J[J0+2] * sigma.r_deexcitation(t,J0) 
        der_sigma_ex_rot      = ( sigma.r_excitation(t,J0)   - sigma.r_excitation(t-dt/beta,J0) ) * beta/dt
        der_sigma_deex_rot    = ( sigma.r_deexcitation(t,J0) - sigma.r_deexcitation(t-dt/beta,J0) ) * beta/dt
        der_total_sigma_term += pN2.x_J[J0] * der_sigma_ex_rot + pN2.x_J[J0+2] * der_sigma_deex_rot
    if switch_vibr != 0:
        for v in range(8):
            for k in range(8):
                if v != k:
                    Ex, Deex = it.vibrational_term(t, f0, v, k, f0D, func)
                    V_Ex    += Ex
                    V_Deex  += Deex
                    total_sigma          += 0.1 * sigma.v_excitation_vk(t, v, k) + 0.1 * sigma.v_deexcitation_vk(t,v,k) 
                    der_sigma_ex_rot      = ( sigma.v_excitation_vk(t,v,k)   - sigma.v_excitation_vk(t-dt/beta,v,k) ) * beta/dt
                    der_sigma_deex_rot    = ( sigma.v_deexcitation_vk(t,v,k) - sigma.v_deexcitation_vk(t-dt/beta,v,k) ) * beta/dt
                    der_total_sigma_term += 0.1 * der_sigma_ex_rot + 0.1 * der_sigma_deex_rot
    total_sigma += sigma.elastic(t)
    der_total_sigma       = (1./total_sigma**2) * (total_sigma - t * der_total_sigma_term)
    A      = t*c.Part_Kn**2/(N_all**2*total_sigma) + (t**2)*c.delta*sigma.elastic(t)
    D      = c.delta*(derivate_term+t**2*sigma.elastic(t))
    C      = c.Part_Kn**2*(der_total_sigma)/(N_all**2)
    res[0] = p # df0/deps
    res[1] = - (1./A)*(p*(D+C) + f0*(derivate_term*c.delta) + (F_ex + F_deex) + switch_vibr*(V_Ex + V_Deex))
    return res

"""Поиск f0"""
y0 = [1, 1.0e-10] #Начальные условия
f0D  = np.zeros(N)+1 #массив значений f0
f_check = np.zeros(N)+5
func = interpolate.UnivariateSpline(electron_energy, f0D, s=0)
epoch = 0
epoch_max = 7
while (epoch < epoch_max):
    f_check = copy.deepcopy(f0D) 
    sol = odeint(search_f0, y0, electron_energy, hmin = dt*1.e-100, hmax=dt*10, mxstep=10**9) #решатель ODE 
    f0D = sol[:,0] #Решение уравнения     
    func = interpolate.UnivariateSpline(electron_energy, f0D, s=0)
#    print(f0D[100]/f_check[100])
    print(epoch, (np.abs((f0D[np.argmax(f_check)]-f_check[np.argmax(f_check)])/f0D[np.argmax(f_check)])))
    epoch += 1
    
for i in range(0,N): 
    f0D[i] *= np.sqrt(electron_energy[i])
#Нормировка интеграла на единицу:
integr = simps(f0D,electron_energy) 
f0D   *= 1./integr
func = interpolate.UnivariateSpline(electron_energy, f0D, s=0) 

#fig = plt.figure(facecolor='white')
plt.plot(electron_energy*c.T/11600., func(electron_energy))
#plt.plot(electron_energy,fD)
plt.legend(('0В, без неупругих', '0В, с неупругими'))
plt.ylabel(r"$f0$")
plt.xlabel(r"$\epsilon$")
plt.grid(True)

working = time.time() - start_time
print("Время работы: %.f минут и %.f секунд" % (int(working/60), working%60))