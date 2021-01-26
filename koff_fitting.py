# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 10:10:52 2020

@author: Allison
"""

# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt

# Liu et. al k_off vs. force graphs
# catch-slip
liu_purple = pd.read_csv("liu_purple.csv")
purple_time = liu_purple['Time']
purple_k_off = 1/purple_time
purple_force = liu_purple['Force']

liu_red = pd.read_csv("liu_red.csv")
red_time = liu_red['Time']
red_k_off = 1/red_time
red_force = liu_red['Force']

liu_blue = pd.read_csv("liu_blue.csv")
blue_time = liu_blue['Time']
blue_k_off = 1/blue_time
blue_force = liu_blue['Force']

# slip
liu_dark_blue = pd.read_csv("liu_dark_blue.csv")
dark_blue_time = liu_dark_blue['Time']
dark_blue_k_off = 1/dark_blue_time
dark_blue_force = liu_dark_blue['Force']

liu_green = pd.read_csv("liu_green.csv")
green_time = liu_green['Time']
green_k_off = 1/green_time
green_force = liu_green['Force']

k_b = 0.0138
temp = 310

def slip(x,y):
    slope, log_k_off_0 = np.polyfit(x,y,1)
    x_B = slope*k_b*temp
    k_off_0 = np.exp(log_k_off_0)
    return x_B, k_off_0

def slip_func(x,x_B,k_off_0):
    return k_off_0*np.exp((x_B*x)/(k_b*temp))

def catch_slip(f,E_21,k_1rup,f_12,k_2rup,x_B):
    exp_1 = np.exp(E_21/(k_b*temp))
    exp_2 = np.exp(f/f_12)
    exp_3 = np.exp((x_B*f)/(k_b*temp))
    return (exp_1*k_1rup + exp_2*k_2rup*exp_3) / (exp_1 + exp_2)

def rsquared_slip(x,y,x_B,k_off_0):
    # yfit = (x_B/(k_b*temp))*x + np.log(k_off_0)
    yfit = k_off_0*np.exp((x_B*x)/(k_b*temp))
    ymean=np.mean(y)
    sstot=sum((y-ymean)**2)
    ssres=sum((y-yfit)**2)
    rs=1-ssres/sstot
    return rs

def rsquared_catch_slip(f,x,y):
    popt,pcov = so.curve_fit(f,x,y,bounds=np.array([0,np.inf]))
    residuals = y - f(x,*popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y-np.mean(y))**2)
    rsq = 1 - (ss_res/ss_tot)
    return rsq

# catch-slip fitting
# fit parameters: [E_21,k_1rup,f_12,k_2rup,x_B]
red1, red2 = so.curve_fit(catch_slip,red_force,red_k_off,bounds=np.array([0,np.inf]))
blue1, blue2 = so.curve_fit(catch_slip,blue_force,blue_k_off,bounds=np.array([0,np.inf]))
purple1, purple2 = so.curve_fit(catch_slip,purple_force,purple_k_off,bounds=np.array([0,np.inf]))
green1, green2 = so.curve_fit(catch_slip,green_force,green_k_off,bounds=np.array([0,np.inf]))
dark_blue1, dark_blue2 = so.curve_fit(catch_slip,dark_blue_force,dark_blue_k_off,bounds=np.array([0,np.inf]))

E_21_list = [red1[0],blue1[0],purple1[0],green1[0],dark_blue1[0]]
k_1rup_list = [red1[1],blue1[1],purple1[1],green1[1],dark_blue1[1]]
f_12_list = [red1[2],blue1[2],purple1[2],green1[2],dark_blue1[2]]
k_2rup_list = [red1[3],blue1[3],purple1[3],green1[3],dark_blue1[3]]
x_B_list_cs = [red1[4],blue1[4],purple1[4],green1[4],dark_blue1[4]]

# catch_slip = cs
cs_rsq_red = rsquared_catch_slip(catch_slip,red_force,red_k_off)
cs_rsq_purple = rsquared_catch_slip(catch_slip,purple_force,purple_k_off)
cs_rsq_blue = rsquared_catch_slip(catch_slip,blue_force,blue_k_off)
cs_rsq_green = rsquared_catch_slip(catch_slip,green_force,green_k_off)
cs_rsq_dark_blue = rsquared_catch_slip(catch_slip,dark_blue_force,dark_blue_k_off)

# slip fitting
x_B1, k_off_01 = slip(red_force,np.log(red_k_off))
rsq_red = rsquared_slip(red_force,np.log(red_k_off),x_B1,k_off_01)
red_s_fit = slip_func(red_force,x_B1,k_off_01)

x_B2, k_off_02 = slip(purple_force,np.log(purple_k_off))
rsq_purple = rsquared_slip(purple_force,np.log(purple_k_off),x_B2,k_off_02)
purple_s_fit = slip_func(purple_force,x_B2,k_off_02)

x_B3, k_off_03 = slip(blue_force,np.log(blue_k_off))
rsq_blue = rsquared_slip(blue_force,np.log(blue_k_off),x_B3,k_off_03)
blue_s_fit = slip_func(blue_force,x_B3,k_off_03)

x_B4, k_off_04 = slip(green_force,np.log(green_k_off))
rsq_green = rsquared_slip(green_force,np.log(green_k_off),x_B4,k_off_04)
green_s_fit = slip_func(green_force,x_B4,k_off_04)

x_B5, k_off_05 = slip(dark_blue_force,np.log(dark_blue_k_off))
rsq_dark_blue = rsquared_slip(dark_blue_force,np.log(dark_blue_k_off),x_B5,k_off_05)
dark_blue_s_fit = slip_func(dark_blue_force,x_B5,k_off_05)

x_B_list_s = [x_B1,x_B2,x_B3,x_B4,x_B5]
k_off_0_list = [k_off_01,k_off_02,k_off_03,k_off_04,k_off_05]
s_fit_list = [red_s_fit,purple_s_fit,blue_s_fit,green_s_fit,dark_blue_s_fit]

peptides_list = ['OVA', 'G4','A2','E1','R4']
graph_colors = list(['r','m','c','g','b'])

force_arr = list([red_force,purple_force,blue_force,green_force,dark_blue_force])
k_off_arr = list([red_k_off,purple_k_off,blue_k_off,green_k_off,dark_blue_k_off])
matrix_arr = list([red1,purple1,blue1,green1,dark_blue1])
# new_force = dark_blue_force[0:4]
# new_k_off = dark_blue_k_off[0:4]
# x_B6, k_off_06 = slip(new_force,np.log(new_k_off))
# new1, new2 = so.curve_fit(catch_slip,new_force,new_k_off,bounds=np.array([0,np.inf]))

# s_x_B_list = [x_B1,x_B2,x_B3,x_B4,x_B5]
# numbers = [1,2,3,4,5]

# comparing x_B values
def compare_x_B(x_B_slip_arr,x_B_cs_arr,min_x_B, max_x_B):
    x_B_cs_final = []
    E_21_final = []
    k_1rup_final = []
    k_2rup_final = []
    f_12_final = []
    
    k_off_0_final = []
    x_B_slip_final = []
    
    model_desc = []
    colors_s = []
    colors_cs = []
    s_peptides = []
    cs_peptides = []
    
    s_forces = []
    cs_forces = []
    s_k_off = []
    cs_k_off = []
    
    s_exp = []
    cs_matrices = []

    if len(x_B_slip_arr) != len(x_B_cs_arr):
        print('Arrays must be equal in length.')
    else:
        for i in range(len(x_B_slip_arr)):
            if (x_B_cs_arr[i] < min_x_B) or (x_B_cs_arr[i] > max_x_B):
                x_B_slip_final.append(x_B_slip_arr[i])
                k_off_0_final.append(k_off_0_list[i])
                model_desc.append('s')
                colors_s.append(graph_colors[i])
                s_peptides.append(peptides_list[i])
                s_forces.append(force_arr[i])
                s_k_off.append(k_off_arr[i])
                s_exp.append(s_fit_list[i])
            else:
                x_B_cs_final.append(x_B_cs_arr[i])
                E_21_final.append(E_21_list[i])
                k_1rup_final.append(k_1rup_list[i])
                k_2rup_final.append(k_2rup_list[i])
                f_12_final.append(f_12_list[i])
                model_desc.append('cs')
                colors_cs.append(graph_colors[i])
                cs_peptides.append(peptides_list[i])
                cs_forces.append(force_arr[i])
                cs_k_off.append(k_off_arr[i])
                cs_matrices.append(matrix_arr[i])
    
    # return all lists (initially empty)
    return x_B_cs_final, E_21_final, k_1rup_final, k_2rup_final, f_12_final, k_off_0_final, x_B_slip_final, model_desc, s_peptides, s_forces, s_k_off, s_exp, cs_peptides, cs_forces, cs_k_off, cs_matrices, colors_s, colors_cs

min_x_B = 10e-3
max_x_B = 10
lists = compare_x_B(x_B_list_s,x_B_list_cs,min_x_B,max_x_B)
x_B_cs, E_21_cs, k_1rup_cs, k_2rup_cs, f_12_final = lists[0:5]
k_off_0_s, x_B_s = lists[5:7]
models = lists[7]
s_peptide_list, s_force_list, s_k_off_list, s_fits = lists[8:12]
cs_peptide_list, cs_force_list, cs_k_off_list, cs_matrices = lists[12:16]
s_color_list, cs_color_list = lists[16:]

# slip model data
slip_data = {}
slip_data['Name'] = s_peptide_list
slip_data['x_B'] = x_B_s
slip_data['k_off_0'] = k_off_0_s
slip_pandas = pd.DataFrame(slip_data)
print(slip_pandas)
print()

# catch-slip model data
cs_data = {}
cs_data['Name'] = cs_peptide_list
cs_data['x_B'] = x_B_cs
cs_data['k_1rup'] = k_1rup_cs
cs_data['k_2rup'] = k_2rup_cs
cs_data['f_12'] = f_12_final
cs_data['E_21'] = E_21_cs
cs_pandas = pd.DataFrame(cs_data)
print(cs_pandas)

# plots
plt.figure(figsize=(8,8))

# slip model plot
plt.subplot(2,1,1)
plt.xlabel('Force')
plt.ylabel('k_off')
for i in range(len(s_force_list)):
    plt.plot(s_force_list[i],s_k_off_list[i],'%so' % s_color_list[i],label='%s' % s_peptide_list[i])
    plt.plot(s_force_list[i],s_fits[i],'k')
plt.legend(loc=2)

# catch-slip model plot
plt.subplot(2,1,2)
plt.xlabel('Force')
plt.ylabel('k_off')
for i in range(len(cs_force_list)):
    plt.plot(cs_force_list[i],cs_k_off_list[i],'%so' % cs_color_list[i],label='%s' % cs_peptide_list[i])
    plt.plot(cs_force_list[i],catch_slip(cs_force_list[i],*cs_matrices[i]),'k')
plt.legend(loc=2)
plt.show()

# # comparing R^2 values
# color_list = ['red','purple','blue','green','dark_blue']
# rsq_slip = [rsq_red,rsq_purple,rsq_blue,rsq_green,rsq_dark_blue]
# rsq_catch_slip = [cs_rsq_red,cs_rsq_purple,cs_rsq_blue,cs_rsq_green,cs_rsq_dark_blue]

# def check_rsq(rsq_slip_arr,catch_slip_rsq_arr):
#     s_rsq_list = []
#     cs_rsq_list = []
#     new_E_21 = []
#     new_k_1rup = []
#     new_k_2rup = []
#     new_f_12 = []
#     new_cs_x_B = []
    
#     new_s_x_B = []
#     new_k_off_0 = []
    
#     model_desc = []
#     s_agonists = []
#     cs_agonists = []
    
#     if len(rsq_slip_arr) != len(catch_slip_rsq_arr):
#         print('Invalid data')
#     else:
#         for i in range(len(rsq_slip_arr)):
#             if rsq_slip_arr[i] < catch_slip_rsq_arr[i]:
#                 cs_rsq_list.append(catch_slip_rsq_arr[i])
#                 new_E_21.append(E_21_list[i])
#                 new_k_1rup.append(k_1rup_list[i])
#                 new_k_2rup.append(k_2rup_list[i])
#                 new_f_12.append(f_12_list[i])
#                 new_cs_x_B.append(x_B_list[i])
#                 model_desc.append('cs')
#                 cs_agonists.append(agonist_list[i])
#             elif rsq_slip_arr[i] > catch_slip_rsq_arr[i]:
#                 s_rsq_list.append(rsq_slip_arr[i])
#                 new_s_x_B.append(s_x_B_list[i])
#                 new_k_off_0.append(k_off_0_list[i])
#                 model_desc.append('s')
#                 s_agonists.append(agonist_list[i])
                
#     return s_rsq_list, cs_rsq_list, new_E_21, new_k_1rup, new_k_2rup, new_f_12, new_cs_x_B, new_s_x_B, new_k_off_0, model_desc, s_agonists, cs_agonists

# many_lists = check_rsq(rsq_slip,rsq_catch_slip)

# plots
# plt.figure(figsize=(10,10))

# plt.subplot(2,1,1)
# plt.plot(purple_force,purple_k_off,'mo',label='G4')
# plt.plot(purple_force,catch_slip(purple_force,*purple1),'k')
# plt.plot(red_force,red_k_off,'ro',label='OVA')
# plt.plot(red_force,catch_slip(red_force,*red1),'r')
# plt.plot(blue_force,blue_k_off,'co',label='A2')
# plt.plot(blue_force,catch_slip(blue_force,*blue1),'b')
# plt.ylabel('Catch-slip k_off')
# plt.legend()

# plt.subplot(2,1,2)
# plt.plot(green_force,green_k_off,'go',label='E1')

# new_force = dark_blue_force[0:4]
# new_k_off = dark_blue_k_off[0:4]

# x_B6, k_off_06 = slip(new_force,np.log(new_k_off))
# new_fit = (x_B6/(k_b*temp))*new_force + np.log(k_off_06)
# rsq_new = rsquared_slip(new_force,np.log(new_k_off),x_B6,k_off_06)
# print(rsq_new)
# print()

# new1, new2 = so.curve_fit(catch_slip,new_force,new_k_off,bounds=np.array([0,np.inf]))
# new_rsq_dark_blue = rsquared_catch_slip(catch_slip,new_force,new_k_off)
# print(new_rsq_dark_blue)

# plt.subplot(3,1,1)
# plt.plot(new_force,new_k_off,'bo')
# plt.plot(new_force,catch_slip(new_force,*new1),'k')

# plt.subplot(3,1,2)
# plt.plot(new_force,np.log(new_k_off),'ro')
# plt.plot(new_force,new_fit,'k')

# plt.subplot(3,1,3)
# plt.plot(new_force,new_k_off,'go')
# plt.plot(new_force,k_off_06*np.exp((x_B6*new_force)/(k_b*temp)),'k')
# plt.show()