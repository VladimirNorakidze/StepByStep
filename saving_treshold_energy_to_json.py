# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 15:41:34 2018

@author: Владимир
"""
import json
import numpy as np

B = 2.49e-4 #Вращательная постоянная для N2 [эВ]
T  = 873. # Температура газа [К]

level_energy = {}
threshold_energy = {} #словарь пороговых энергий
threshold_energy["rot"] = []
g = np.ndarray(200) #стат веса 
J1 = 0
level_energy["rot"] = []

for J in range(300):
    level_energy["rot"].append(B*J*(J+1)*11600./T) #Безразмерная энергия уровня J

for J in range(250):
    threshold_energy["rot"].append(level_energy["rot"][J+2]-level_energy["rot"][J])

"""Пороговые энергии колебательных переходов"""
threshold_energy["vibr"] = [[1.8, 1.8, 1.8, 1.912, 2.079, 2.1, 2.299, 2.397, 2.594],
                                  [0, 1.51, 1.51, 1.622, 1.789, 1.81, 2.009, 2.107, 2.304],
                                  [0, 0, 1.21, 1.322, 1.489, 1.51, 1.709, 1.807, 2.004], 
                                  [0, 0, 0, 1.32, 1.199, 1.22, 1.419, 1.517, 1.714],
                                  [0, 0, 0, 0, 0.909, 0.93, 1.129, 1.227, 1.424], 
                                  [0, 0, 0, 0, 0, 0.63, 0.829, 0.927, 1.124], 
                                  [0, 0, 0, 0, 0, 0, 0.539, 0.637, 0.834], 
                                  [0, 0, 0, 0, 0, 0, 0, 0.337, 0.534], 
                                  [0, 0, 0, 0, 0, 0, 0, 0, 0.244]]
for v in range(len(threshold_energy["vibr"])):
       for i in range(len(threshold_energy["vibr"][v])):
           threshold_energy["vibr"][v][i] *=  11600./T 
           
#threshold_energy["vibrational"][0][0] = 1.75
#threshold_energy["vibrational"][0][1] = 1.77
#threshold_energy["vibrational"][0][2] = 1.82
#threshold_energy["vibrational"][0][3] = 1.86 
#threshold_energy["vibrational"][0][4] = 2.02
#threshold_energy["vibrational"][0][5] = 2.08
#threshold_energy["vibrational"][0][6] = 2.27
#threshold_energy["vibrational"][0][7] = 2.47
#threshold_energy["vibrational"][0][8] = 2.56
#threshold_energy["vibrational"][1][0] = 1.46
#threshold_energy["vibrational"][1][1] = 1.48
#threshold_energy["vibrational"][1][2] = 1.53
#threshold_energy["vibrational"][1][3] = 1.57
#threshold_energy["vibrational"][1][4] = 1.73
#threshold_energy["vibrational"][1][5] = 1.79
#threshold_energy["vibrational"][1][6] = 1.98
#threshold_energy["vibrational"][1][7] = 2.18
#threshold_energy["vibrational"][1][8] = 2.27
#threshold_energy["vibrational"][2][0] = 1.18
#threshold_energy["vibrational"][2][1] = 1.2
#threshold_energy["vibrational"][2][2] = 1.25
#threshold_energy["vibrational"][2][3] = 1.29
#threshold_energy["vibrational"][2][4] = 1.45
#threshold_energy["vibrational"][2][5] = 1.51
#threshold_energy["vibrational"][2][6] = 1.7
#threshold_energy["vibrational"][2][7] = 1.9
#threshold_energy["vibrational"][2][8] = 1.99
#threshold_energy["vibrational"][3][0] = 0.89
#threshold_energy["vibrational"][3][1] = 0.91
#threshold_energy["vibrational"][3][2] = 0.96
#threshold_energy["vibrational"][3][3] = 1.
#threshold_energy["vibrational"][3][4] = 1.16
#threshold_energy["vibrational"][3][5] = 1.22
#threshold_energy["vibrational"][3][6] = 1.41
#threshold_energy["vibrational"][3][7] = 1.61
#threshold_energy["vibrational"][3][8] = 1.7
#threshold_energy["vibrational"][4][0] = 0.62
#threshold_energy["vibrational"][4][1] = 0.64
#threshold_energy["vibrational"][4][2] = 0.69
#threshold_energy["vibrational"][4][3] = 0.73
#threshold_energy["vibrational"][4][4] = 0.89
#threshold_energy["vibrational"][4][5] = 0.95
#threshold_energy["vibrational"][4][6] = 1.14
#threshold_energy["vibrational"][4][7] = 1.34
#threshold_energy["vibrational"][4][8] = 1.43
#threshold_energy["vibrational"][5][0] = 0.34
#threshold_energy["vibrational"][5][1] = 0.36
#threshold_energy["vibrational"][5][2] = 0.41
#threshold_energy["vibrational"][5][3] = 0.45
#threshold_energy["vibrational"][5][4] = 0.61
#threshold_energy["vibrational"][5][5] = 0.67
#threshold_energy["vibrational"][5][6] = 0.86
#threshold_energy["vibrational"][5][7] = 1.06
#threshold_energy["vibrational"][5][8] = 1.15
#threshold_energy["vibrational"][6][0] = 0.07
#threshold_energy["vibrational"][6][1] = 0.09
#threshold_energy["vibrational"][6][2] = 0.14
#threshold_energy["vibrational"][6][3] = 0.18
#threshold_energy["vibrational"][6][4] = 0.34
#threshold_energy["vibrational"][6][5] = 0.4
#threshold_energy["vibrational"][6][6] = 0.59
#threshold_energy["vibrational"][6][7] = 0.79
#threshold_energy["vibrational"][6][8] = 0.88
#threshold_energy["vibrational"][7][4] = 0.07
#threshold_energy["vibrational"][7][5] = 0.13
#threshold_energy["vibrational"][7][6] = 0.32
#threshold_energy["vibrational"][7][7] = 0.52
#threshold_energy["vibrational"][7][8] = 0.61
#threshold_energy["vibrational"][8][6] = 0.06
#threshold_energy["vibrational"][8][7] = 0.26
#threshold_energy["vibrational"][8][8] = 0.35

with open("threshold_energy.txt","w",encoding="utf-8") as file:
    json.dump(threshold_energy,file)
    
file.close()