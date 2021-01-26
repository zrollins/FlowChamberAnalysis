# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 08:40:36 2020

@author: Allison
"""
# test data: Brunk et. al (2 graphs)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# u_f = y*gamma_w*(1-(5/16)*(a/y)**3) = U_hd

# variables
a = 5  # cell radius (um)
y = 5.025 # distance from wall to edge of cell (um)
L = 20 # bond length (nm)
w = 800 # chamber width (um)
b = 50 # chamber half height (um)
mu = 6.92e-3 # (dyne/cm^2)

def k_plus_func(U_hd,U_cell,k_off,m_l):
    ratio = U_hd/U_cell
    num = ratio*k_off*(1-(1/ratio))
    return num/m_l

def U_hd_func(y,gamma_w,a):
    return y*gamma_w*(1-(5/16)*(a/y)**3)
    
def tau_func(mu,Q,w,b):
    return (3*mu*Q)/(2*w*b**2)

def gamma_func(a,mu,tau):
    return 0.9440*(4*np.pi*mu*a**3)*(tau/mu)

# def U_cell_func(k_off_arr,U_hd,k_plus):
#     return U_hd*(1-(k_plus*3312)/(k_plus*3312 + k_off_arr))

data_files = []
k_off_files = []
Q_list = []
m_l_list = []
tau_list = []
gamma_list = []
U_hd_list = []

# parameters from Brunk
L = 6.6 # cm
w = 0.1 # cm
b_double = 175 # um
Q = 128 # mL/um
mu = 8.9e-3 # dyne*s/cm^2 (same as water at room temp)
a = 5 # um
m_l = 3312
# use my parameters for y
# tau = gamma*mu

# constant parameters
# mu = float(input('Enter viscosity (dyne*s/cm^2): '))
# a = float(input('Enter cell radius (um): '))
# y = float(input('Enter distance from flow chamber to edge of cell (um): '))
# L = float(input('Enter bond length (nm): '))
# w = float(input('Enter flow chamber width (um): '))
# b_double = float(input('Enter flow chamber height (um): '))
# b = b_double/2

# user's input for now:
# site density for that file
# parameters to calculate U_f = U_hd
# file with U_cell vs. tau data
# Q to calculate tau

# imports:
# k_off (from part 1); from part_1 import *

while True:
    run = str(input('Enter \"y\" to input files and flow rate or \"n\" to quit: '))
    
    if run.lower() == 'n':
        break
    
    elif run.lower() == 'y':
        
        # flow rate
        Q = float(input('Enter flow rate: '))
        Q_list.append(Q)
        
        # site density
        m_l = float(input('Enter site density: '))
        m_l_list.append(m_l)
        
        # k_off
        # use fitting parameters to calculate k_off (slip/catch-slip)
        # use papers 14, 24 (Cheung)
        # validate using previous info (papers in Cheung)
        # use low force (< 0.5 dyn/cm^2)
        # 0.5: theoretical diverging point for secondary k+
        # plot k+ vs tau
        # plot both k+ curves on same graph, compare them (should be the same?)
        # show that both curves are qualitatively similar
        # compare eqn 8 and 10 (Cheung) for validation
        # include R^2 just to show that data fits model
        k_off_data = input('Enter names of files separated by commas on one line with bond lifetimes and k_off data: ')
        result = [x.strip() for x in k_off_data.split(',')]
        for i in range(len(result)):
            if '.csv' not in result[i]:
                result[i] += '.csv'
        
        try:
            for i in range(len(result)):
                with open(k_off_data) as file_open:
                    file = file_open.read()
                    k_off_files.append(result[i])
        except FileNotFoundError:
            print('Invalid file name.')
        
        # gamma
        tau = tau_func(mu,Q,w,b)
        tau_list.append(tau)
        
        gamma = gamma_func(a,mu,tau)
        gamma_list.append(gamma)
        
        U_hd = U_hd_func(y, gamma, a)
        U_hd_list.append(U_hd)
        
        # tracks data
        file_input = input('Enter names of files separated by commas on one line with U_cell vs. shear stress data: ')
        result = [x.strip() for x in file_input.split(',')]
        
        for i in range(len(result)):
            if '.csv' not in result[i]:
                result[i] += '.csv'
        
        try:
            for i in range(len(result)):
                with open(file_input) as file_open:
                    file = file_open.read()
                    data_files.append(result[i])
        except FileNotFoundError:
            print('Invalid file name.')
            
            
    else:
        print('Please enter \"y\" or \"n\".')

panda_bois = []
tau_data = []
U_cell_data = []
k_off_arr = []

# U_cell, tau from input files
for i in range(len(data_files)):
    panda_boi = pd.read_csv(data_files[i])
    panda_bois.append(panda_boi)
    
    tau_boi = panda_boi['Tau']
    tau_data.append(tau_boi)
    
    U_cell_boi = panda_boi['U_cell']
    U_cell_data.append(U_cell_boi)

U_cell_arr = np.zeros(len(result))
tau_arr = np.zeros(len(result))

for i in range(len(tau_data)):
    for j in range(len(result)):
        U_cell_arr[j] = U_cell_data[i][j]
        tau_arr[j] = tau_data[i][j]

# k_off from input files
for i in range(len(k_off_files)):
    panda_boi = pd.read_excel(k_off_files[i])
    k_off_boi = panda_boi['k_off']
    k_off_arr.append(k_off_boi)
    
U_hd_arr = np.array(U_hd_list)

k_plus_lists = []

for i in range(len(Q_list)):
    k_plus_list = []
    for j in range(len(U_cell_arr)):
        k_plus_list.append(k_plus_func(U_hd_arr[i], U_cell_arr[j], k_off_arr[i][j], m_l_list[i]))   
    k_plus_lists.append(k_plus_list)
    # plot should look like line in Fig 3B (Cheung)

# # U_cell_list = []
# # for i in range(len(tracks_data_files)):
# #     tracks_data = pd.read_csv(tracks_data_files[i])
# #     tracks_speed = tracks_data['TRACK_MEAN_SPEED']
# #     U_cell_list.append(np.mean(tracks_speed))

# U_cell_arr = U_cell_func(U_hd,k_plus,k_off,50)
# plt.plot(tau,U_cell_arr,'bo')
# plt.show()