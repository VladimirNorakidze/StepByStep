# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 13:50:24 2018

@author: Владимир
"""

import sigma 
import json

"""Загрузка пороговых энергий неупругих переходов молекулы азота"""
with open("threshold_energy.txt","r") as file:
    threshold_energy = json.load(file)
    
