# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:31:15 2020

@author: Allison
"""
import numpy as np
from scipy import optimize

m_r = int(input('Enter receptor ligand density (1/um^2): ')) # m_r = 72
Nb = input('Enter Nb: ')
Nb_list = []
Nb_str = [val.strip() for val in Nb.split(',')]
for i in range(len(Nb_str)):
    Nb_list.append(int(Nb_str[i]))
m_l = str(input('Enter site densities (sites/um^2) separated by commas and without new lines: '))
m_l_list = []
m_l_str = [val.strip() for val in m_l.split(',')]
for i in range(len(m_l_str)):
    m_l_list.append(int(m_l_str[i]))
    
mrml = m_r*m_l_list

def Nb_func(mrml,NT,AcKa):
    dem = AcKa*mrml
    return NT/(1+(1/dem))

parameters, matrix = optimize.curve_fit(Nb_func,mrml,Nb_list,bounds=np.array([0,np.inf]))
NT = parameters[0]