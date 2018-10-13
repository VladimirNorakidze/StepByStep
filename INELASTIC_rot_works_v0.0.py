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

import sigma
import rotational_term as rt
#import main_function as mf
import population_N2 as pN2

from numba import vectorize

vectorize(['float32(float32, float32)'], target='cuda')

"""
Constants
"""
e  = 1.6*10**(-19) #Кл
m  = 9.1*10**(-31) # кг - масса электрона
M0 = 1.67*10**(-27) #кг - масса нуклона
k  = 1.38e-23 #Постоянная Больцмана
c = 3.0e8 #Сскорость света (м/c)
Q      = -1.13 #Квадрапольный момент молекулы N2 (в единиицах e*a0^2)
alpha   = 3.13 #поляризуемость (в единицах a0^3)
Planck_const = 1.0545718e-34 #Дж*с
a0     = 137.*Planck_const/(m*c) #Радиус Бора (м)
B = 2.49e-4 #Вращательная постоянная для N2 [эВ]
diss_energy_N2 = 9.73 # [эВ]
"""
Parameters
"""
A  = 28 #Массовое число
U  = 630. #Напряжение [В]
r  = 0.9 #Расстояние между электродами [м]
E  = U/r #[В/м]
M  = A*M0 #масса ядра [кг]
delta=2.*m/M
p  = 1.5*133 #Давление Torr*(Па/Torr) (итоговая размерность [Па])
T  = 873. # Температура газа [К]
n0 = 1.0e16 #Плотность электронов [м^(-3)]
sigma0 = np.pi*a0**2 # Площадь боровской орбиты [м^2]
Part_Kn = np.sqrt(1./3)*e*E/( (k*T) * sigma0) #Kn=Part_Kn/(N_all*sigma*(корень из eps))
beta = 1.0e8

N  = 13000 #число точек
N_for_dif_case = 1.3e4
electron_energy = np.linspace( 1.0e-4, N_for_dif_case, N) # значения энергии

with open("constants.py", 'w') as f:
    f.write('import numpy as np \n\n\
e = %.20f\n\
m = %.32f\n\
M0 = %.29f\n\
k = %.25f\n\
c = %f\n\
Q = %f\n\
alpha = %f\n\
Planck_const = %.40f\n\
a0 = 137.*Planck_const/(m*c)\n\
B = %f\n\
diss_energy_N2 = %f\n\n\
A = %f\n\
U = %f\n\
r = %f\n\
E = U/r\n\
M = A*M0\n\
delta = 2.*m/M\n\
p = %f\n\
T = %f\n\
n0 = %f\n\
sigma0 = np.pi*a0**2\n\
Part_Kn = np.sqrt(1./3)*e*E/( (k*T) * sigma0)\n\
beta = %f' %(e,m,M0,k,c,Q,alpha,Planck_const,B,diss_energy_N2,A,U,r,p,T,n0,beta))


"""Расчет числа молекул. Они подчиняются закону распределения Максвелла"""
N_all  = p/(k*T) #м^-3

"""Поиск заселенностей основного J0=0 и возбужденных уровней"""
g, x_J, eps_J, eps_th_ex = pN2.population_N2(B, diss_energy_N2)
"""Сечения возбуждения колебательных степеней"""
eps  = np.array([0, 0.5, 1.0, 1.5, 1.98, 2.1, 2.46, 2.605, 3.0, 5.0, 7.5, 10.0, 
                 15.0, 18.0, 20.0, 22.5, 25.0, 30.0, 1e17])*(11600./T)
    
sig  = np.array([0, 0.005, 0.009, 0.089, 4.56, 1.97, 1.65, 4.4, 1.37, 0.08, 0.031, 
                 0.015, 0.039, 0.076, 0.195, 0.126, 0.082, 0.027 ,0])/(sigma0*1.0e20)

"""Функция сечения возбуждения колебательных степеней для перехода 0->1:"""
sigma_vibr_ex = interpolate.UnivariateSpline(eps, sig, s=0)

#def sigma_vibr_deex(energy,v=0,k=1): #Функция для перехода k->v
#    x_kv = ( energy + 0.289 )/energy #0.289 - Энергия перехода 0 -> 1
#    sigma = sigma_vibr_ex(energy + 0.289) * x_kv
#    return sigma/sigma0
#
#
#def vibrational_term(eps, f0):
#    if (f0D[10] == 1):
#        F_ex   = -1   * f0 * (eps*sigma_vibr_ex(eps)   - (eps+0.289)*sigma_vibr_ex(eps+0.289))
#        F_deex = -1 * f0 * (eps*sigma_vibr_deex(eps) - (eps-0.289)*sigma_vibr_deex(eps-0.289))
#    else:
#        f1 = func(eps+0.289)
#        if eps-0.289 <= 0:
#            f2 = 0
#        else:
#            f2 = func(eps-0.289)
#        F_ex   = -1   * (eps*sigma_vibr_ex(eps)*func(eps)   - (eps+0.289)*sigma_vibr_ex(eps+0.289)*f1)
#        F_deex = -1 * (eps*sigma_vibr_deex(eps)*func(eps) - (eps-0.289)*sigma_vibr_deex(eps-0.289)*f2)
#    return [F_ex, F_deex] 
#count = 0
#NN = 0

"""Функция для решения дифференцияального уравнения относительно f0:"""    
def search_f0(y,t):
    res   = np.ndarray((2))
    f0, p = y
    total_sigma = 0
    derivate_term = t*(2*sigma.elastic(t) + t*derivative(sigma.elastic,t))
    der_total_sigma_term = derivative(sigma.elastic,t)
    F_ex   = 0
    F_deex = 0
#    J0     = 0
    for J0 in range(pN2.N_rot_trans(B, diss_energy_N2)):
        Ex, Deex = rt.rotational_term(t, J0, x_J, f0, f0D, eps_th_ex[J0], func)
        F_ex   += Ex
        F_deex += Deex
        total_sigma          += x_J[J0] * sigma.r_excitation(t,J0, eps_th_ex[J0]) + x_J[J0+2] * sigma.r_deexcitation(t,J0, eps_th_ex[J0]) 
        der_sigma_ex_rot      = ( sigma.r_excitation(t,J0, eps_th_ex[J0])   - sigma.r_excitation(t-dt/beta,J0, eps_th_ex[J0]) ) * beta/dt
        der_sigma_deex_rot    = ( sigma.r_deexcitation(t,J0, eps_th_ex[J0]) - sigma.r_deexcitation(t-dt/beta,J0, eps_th_ex[J0]) ) * beta/dt
        der_total_sigma_term += x_J[J0] * der_sigma_ex_rot + x_J[J0+2] * der_sigma_deex_rot
    V_Ex, V_Deex = 0, 0
    total_sigma += sigma.elastic(t)
    der_total_sigma       = (1./total_sigma**2) * (total_sigma - t * der_total_sigma_term)
    A      = t*Part_Kn**2/(N_all**2*total_sigma) + (t**2)*delta*sigma.elastic(t)
    D      = delta*(derivate_term+t**2*sigma.elastic(t))
    C      = Part_Kn**2*(der_total_sigma)/(N_all**2)
    res[0] = p # df0/deps
    res[1] = - (1./A)*(p*(D+C) + f0*(derivate_term*delta) + (F_ex + F_deex) + (V_Ex + V_Deex))
    return res

"""Поиск f0"""
dt   = electron_energy[10] - electron_energy[9] #Шаг по энергии 
y0 = [1, 1.0e-10] #Начальные условия
f0D  = np.zeros(N)+1 #массив значений f0
f_check = np.zeros(N)+5
func = interpolate.UnivariateSpline(electron_energy, f0D, s=0)
z = 0
while (z < 1):
    f_check = copy.deepcopy(f0D) 
    sol = odeint(search_f0, y0, electron_energy, hmin = dt*1.e-100, hmax=dt*10, mxstep=10**9) #решатель ODE 
    f0D = sol[:,0] #Решение уравнения     
    func = interpolate.UnivariateSpline(electron_energy, f0D, s=0)
#    print(f0D[100]/f_check[100])
    print(z, (np.abs((f0D[100]-f_check[100])/f0D[100])))
    z+=1
for i in range(0,N): 
    f0D[i] *= np.sqrt(electron_energy[i])
#Нормировка интеграла на единицу:
integr = simps(f0D,electron_energy) 
f0D   *= 1./integr
func = interpolate.UnivariateSpline(electron_energy, f0D, s=0) 

#fig = plt.figure(facecolor='white')
plt.plot(electron_energy*T/11600., func(electron_energy))
#plt.plot(electron_energy,fD)
plt.legend(('0В, без неупругих', '0В, с неупругими'))
plt.ylabel(r"$f0$")
plt.xlabel(r"$\epsilon$")
plt.grid(True)

