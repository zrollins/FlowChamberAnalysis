# -*- coding: utf-8 -*-
"""
Created on Mon May 31 19:11:53 2021

@author: Allison
"""

import pandas as pd
import numpy as np
from scipy import stats, optimize, special, interpolate
import scipy.stats as stats
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.special as sp
import xlsxwriter
import operator


import os
os.chdir('C:/Users/Allison/Downloads/lab')

bob = pd.read_csv('control_spots2.csv')

# %%
# # variables (should be inputted by user)
# a = 5  # cell radius (um)
# y = 5.025 # distance from wall to edge of cell (um)
# gamma_w = 1 # wall shear rate 
# L = 20 # bond length (nm)
# w = 800 # chamber width (um)
# b = 50 # chamber half height (um)
# mu = 6.92e-3 # (dyne/cm^2)

# empty lists to store data
parameters = []

# experimental parameters
# assume these remain constant
mu = float(input('Enter viscosity (dyne-s/cm^2): ')) # check viscosity units
a = float(input('Enter cell radius (microns): '))
d = float(input('Enter critical distance (microns): '))
L = float(input('Enter receptor-ligand bond length (nm): '))
b = float(input('Enter flow chamber height (microns): '))
b /= 2
w = float(input('Enter flow chamber width (microns): '))
parameters.append([mu, a, b, L, w, d])

y = a+d

CCD_FPS = int(input('Enter CCD FPS: '))

# inter/extrapolation
da = np.array([0, 10e-8, 10e-7, 10e-6, 
               10e-5, 10e-4, 10e-3, 0.003202])
speed_constants = np.array([0.5676, 0.5556, 0.5539, 0.5518, 
                           0.5488, 0.5445, 0.5375, 0.5315])
fit_func = interpolate.interp1d(da, speed_constants,
                                fill_value='extrapolate')

speed_const = fit_func(d/a)

# %% parameters with unit conversions
# mu = float(input('Enter viscosity (dyne-s/cm^2): ')) * 1e14
# a = float(input('Enter cell radius (microns): ')) * 1e6
# d = float(input('Enter critical distance (microns): ')) * 1e6
# L = float(input('Enter receptor-ligand bond length (nm): ')) * 1e3
# b = float(input('Enter flow chamber height (microns): ')) * 10e6
# b /= 2
# w = float(input('Enter flow chamber width (microns): ')) * 10e6
# parameters.append([mu, a, b, L, w, d])

# y = a+d

# CCD_FPS = int(input('Enter CCD FPS: '))

# # inter/extrapolation
# da = np.array([0, 10e-8, 10e-7, 10e-6, 
#                10e-5, 10e-4, 10e-3, 0.003202])
# speed_constants = np.array([0.5676, 0.5556, 0.5539, 0.5518, 
#                            0.5488, 0.5445, 0.5375, 0.5315])
# fit_func = interpolate.interp1d(da, speed_constants,
#                                 fill_value='extrapolate')

# speed_const = fit_func(d/a)

# %% loop that inputs user data
track_data = []
spots_data = []
site_densities = []
shear_rates = []
forces = []
t_min_vals = []

while True:
    run = str(input('Enter \"y\" to input data or \"n\" to quit: '))
    
    if run.lower() == 'n':
        break
    
    elif run.lower() == 'y':
        
        # site density
        m_l = int(input('Enter site density (sites/um^2): '))
        site_densities.append(m_l)
        
        # flow rate
        Q = input('For site density = %d, enter flow rates, separated by commas and without new lines: ' % m_l)
        Q_str = [val.strip() for val in Q.split(',')]
        
        Q_arr = np.zeros(len(Q_str))
        for i in range(len(Q_str)):
            Q_arr[i] = int(Q_str[i])
            
        # tether force
        f = Q_arr * np.sqrt(a/(2*L)) * (1.7005*9*np.pi*mu*a**2 + 0.9440*6*np.pi*mu*a**2) / (w*b**2)
        forces.append(list(f))
        
        # shear stress
        tau = (3*mu*Q_arr) / (2*w*b**2)
        shear_rate = tau / mu
        shear_rates.append(shear_rate)
        
        # Trackmate files
        track_data_sublist = []
        spots_data_sublist = []
        t_min_sublist = []
        
        for i in range(len(Q_arr)):
            track_file_name = input('For flow rate = %.2f and site density = %d, enter name of "Track statistics" file (s) from Trackmate: ' % (Q_arr[i], m_l))
            track_file_list = [val.strip() for val in track_file_name.split(',')]
            track_file_subsub = []
            for j in range(len(track_file_list)):
                if '.csv' not in track_file_list[j]:
                    track_file_list[j] += '.csv'
                
                try:
                    with open(track_file_list[j]) as file_open:
                        file = file_open.read()
                        track_file_subsub.append(track_file_list[j])
                
                except FileNotFoundError:
                    print('Invalid file name.')
                    
            track_data_sublist.append(track_file_subsub)
                
            spots_file_name = input('For flow rate = %.2f and site density = %d, enter name of "Spots in track statistics" file (s) from Trackmate: ' % (Q_arr[i], m_l))
            spots_file_list = [val.strip() for val in spots_file_name.split(',')]
            spots_file_subsub = []
            for j in range(len(spots_file_list)):
                if '.csv' not in spots_file_list[j]:
                    spots_file_list[j] += '.csv'
                    
                try:
                    with open(spots_file_list[j]) as file_open:
                        file = file_open.read()
                        spots_file_subsub.append(spots_file_list[j])
                
                except FileNotFoundError:
                    print('Invalid file name.')
                    
            spots_data_sublist.append(spots_file_subsub)
            
            t_min_input = input('If applicable, for flow rate = %.2f and site density = %d, enter non-specific binding time(s). Otherwise, enter \"n\" for all trials: ' % (Q_arr[i], m_l))
            t_min_str = [val.strip() for val in t_min_input.split(',')]
            t_min_subsub = []
            
            for i in range(len(t_min_str)):
                if t_min_str[i] == 'n':
                    t_min_subsub.append(0.2)
                
                else:
                    t_min_float = float(t_min_str[i])
                    t_min_subsub.append(t_min_float)
                    
            t_min_sublist.append(t_min_subsub)
        
        track_data.append(track_data_sublist)
        spots_data.append(spots_data_sublist)
        t_min_vals.append(t_min_sublist)
                 
    else:
        print('Please enter \"y\" or \"n\".')
        
# %% calculating k_off
koff_avg_vals = []
koff_error_vals = []
Nb_vals = []
NbNT_error_vals = []
NbNT_vals = []
u_f_vals = []

total_track_durations = []
total_cell_bound_time = []
trackIDs_unique = []

# for testing only!
shear_rates = [[1,1,1]]

for m in range(len(track_data)):
    
    koff_avg_vals_sublist = []
    koff_error_vals_sublist = []
    Nb_vals_sublist = []
    NbNT_error_vals_sublist = []
    NbNT_vals_sublist = []
    u_f_sublist = []
    
    for n in range(len(track_data[m])):
        
        koff_avg_vals_subsub = []
        koff_error_vals_subsub = []
        Nb_vals_subsub = []
        NbNT_vals_subsub = []
        
        # cell velocity filtering
        u_f = y*shear_rates[m][n]*(1-(5/16)*(a/y)**3)
        u_f_sublist.append(u_f)
        
        # u_f_actual = y*shear_rates[m][n]*(1-(5/16)*(a/y)**3)
        # u_f_actual_sublist.append(u_f_actual)
        
        for p in range(len(track_data[m][n])):
            tracks_raw_data = pd.read_csv(track_data[m][n][p])
            spots_raw_data = pd.read_csv(spots_data[m][n][p])
            
            filtered_speeds = tracks_raw_data[tracks_raw_data['TRACK_MEAN_SPEED'] < np.absolute(u_f)]
            filtered_tracks_list = list(filtered_speeds['TRACK_ID'])
            
            # only collect track IDs from spots stats file for certain velocities
            better_tracks = []
            
            # obtaining track ID's present in both spots and track stats spreadsheets
            trackID = spots_raw_data['TRACK_ID']
            particleID = spots_raw_data['ID']
            x_pos = spots_raw_data['POSITION_X']
            y_pos = spots_raw_data['POSITION_Y']
            frame = spots_raw_data['FRAME']
            for i in range(len(filtered_tracks_list)): 
                for j in range(len(trackID)):
                    if trackID[j] == filtered_tracks_list[i]:
                        if j != 0:
                            if trackID[j-1] != trackID[j]:
                                better_tracks.append(trackID[j])
                        else:
                            better_tracks.append(trackID[j])
                            
            # new lists to categorize data after velocity filtering
            particleID_new = []
            trackID_new = []
            x_new = []
            y_new = []
            frame_new = [] 
            
            # adding better_tracks corresponding data to empty lists 
            for i in range(len(better_tracks)):
                for j in range(len(trackID)):
                    if trackID[j] == better_tracks[i]:
                        particleID_new.append(particleID[j])
                        trackID_new.append(trackID[j])
                        x_new.append(x_pos[j])
                        y_new.append(y_pos[j])
                        frame_new.append(frame[j])
            
            # r refers to meeting criteria
            # find which particles meet stopping criteria
            r_pos_x = []
            r_pos_y = []
            r_trackID = []
            r_particleID = []
            r_frame = []
            
            i = 0
            i_max = len(trackID_new)
            j = 0
            
            # calculate displacement of a given cell
            def calc_disp(x0,x,y0,y):
                return np.sqrt((x-x0)**2+(y-y0)**2)
            
            # filter using stopping criteria
            tmin_frames = t_min_vals[m][n][p] * CCD_FPS
            while i < i_max-1:
                disp1 = calc_disp(x_new[i+1],x_new[j],y_new[i+1],y_new[j])
                if disp1 <= 1:
                    i += 1
                    disp2 = calc_disp(x_new[i],x_new[j],y_new[i],y_new[j])
                    if i-j > tmin_frames:
                        r_particleID.append(particleID_new[i])
                        r_trackID.append(trackID_new[i])
                        r_pos_x.append(x_new[i])
                        r_pos_y.append(y_new[i])
                        r_frame.append(frame_new[i])
                else:
                    i += 1
                    j = i-1
            
            # time conversion: (# of frames) -> seconds
            # tc = time conversion
            tc_particleID = np.array(r_particleID)
            tc_trackID = np.array(r_trackID)
            tc_frame = np.array(r_frame)
            
            # initial parameters
            i = 1
            j = 0
            k = 0
            t_total = []
            tc_particleID_new = []
            tc_trackID_new = []
            t_tot = 0
            
            # doing time conversion
            while i < len(tc_trackID):
                if tc_trackID[i] == tc_trackID[j]:
                    if tc_frame[i]-tc_frame[k] == 1:
                        t_tot += (tc_frame[i] - tc_frame[k])
                        
                        if i == len(tc_trackID)-1:
                            t_total.append((t_tot + tmin_frames + 1) / CCD_FPS)
                            tc_trackID_new.append(tc_trackID[k])
                        
                    else:
                        if i == len(tc_trackID)-1:
                            t_total.append((t_tot + tmin_frames + 1) / CCD_FPS)
                            tc_trackID_new.append(tc_trackID[k])
                            
                            t_total.append((tmin_frames + 1) / CCD_FPS)
                            tc_trackID_new.append(tc_trackID[i])
                            
                        else:
                            t_total.append((t_tot + tmin_frames + 1) / CCD_FPS)
                            tc_trackID_new.append(tc_trackID[k])
                            
                        t_tot = 0
                        j = i
                        
                        
                    i += 1
                    k += 1
                
                else:
                    t_total.append((t_tot + tmin_frames + 1) / CCD_FPS)
                    tc_trackID_new.append(tc_trackID[k])
                    
                    if i == len(tc_trackID) - 1:
                        t_total.append((tmin_frames + 1) / CCD_FPS)
                        tc_trackID_new.append(tc_trackID[i])
                        
                    t_tot = 0
                    j = i
                    i += 1
                    k += 1
                    
            # getting list with only unique track IDs
            i = 1
            j = 0
            k = 0
            t_total_unique = []
            tc_trackID_unique = []
            t_tot = 0
            
            while i < len(tc_trackID_new):
                if tc_trackID_new[i] != tc_trackID_new[j]:
                    tc_trackID_unique.append(tc_trackID_new[j])
                    t_tot += t_total[k]
                    t_total_unique.append(t_tot)
                    
                    if i == len(tc_trackID_new) - 1:
                        tc_trackID_unique.append(tc_trackID_new[i])
                        t_total_unique.append(t_total[i])
                        
                    j = i
                    i += 1
                    k += 1
                    t_tot = 0
                    
                else:
                    t_tot += t_total[k]
                    
                    if i == len(tc_trackID_new) - 1:
                        t_tot += t_total[i]
                        tc_trackID_unique.append(tc_trackID_new[i])
                        t_total_unique.append(t_tot)
                        
                    i += 1
                    k += 1
            
            # computing k_off = 1/lifetimes
            k_off_sublist = []
            for p in range(len(t_total_unique)):
                if t_total_unique[p] != 0:
                    k_off_sublist.append(1/t_total_unique[p])
                else:
                    k_off_sublist.append(0) 
                    # if lifetime = 0, just append 0 as a placeholder
                    # there probably shouldn't be lifetimes of 0 anyway?
                    # maybe add something in stopping criteria loop to prevent lifetime = 0
            
            # getting entire k_off arrays into list
            avg_k_off = np.mean(k_off_sublist)
            # std_error = stats.sem(k_off_sublist)
            Nb = len(tc_trackID_unique)
            NT = len(trackID_new)
            
            koff_avg_vals_subsub.append(avg_k_off)
            # koff_error_vals_subsub.append(std_error)
            Nb_vals_subsub.append(Nb)
            NbNT_vals_subsub.append(Nb/NT)
        
        koff_avg_new = np.mean(koff_avg_vals_subsub)
        koff_avg_error = stats.sem(koff_avg_vals_subsub)
        Nb_avg_new = np.mean(Nb_vals_subsub)
        NbNT_avg_new = np.mean(NbNT_vals_subsub)
        NbNT_avg_error = stats.sem(NbNT_vals_subsub)
        
        koff_avg_vals_sublist.append(koff_avg_new)
        koff_error_vals_sublist.append(koff_avg_error)
        Nb_vals_sublist.append(Nb_avg_new)
        NbNT_error_vals_sublist.append(NbNT_avg_error)
        NbNT_vals_sublist.append(NbNT_avg_new)
        
    koff_avg_vals.append(koff_avg_vals_sublist)
    koff_error_vals.append(koff_error_vals_sublist)
    Nb_vals.append(Nb_vals_sublist)
    NbNT_error_vals.append(NbNT_error_vals_sublist)
    NbNT_vals.append(NbNT_vals_sublist)
    
# %% k_off fitting
k_b = 0.0138
temp = 310 # Kelvin

# for testing only
forces = [[4, 6, 8]]

# slip model fitting parameters
def slip(x,y):
    slope, log_k_off_0 = np.polyfit(x,y,1)
    x_B = slope*k_b*temp
    k_off_0 = np.exp(log_k_off_0)
    return x_B, k_off_0

# slip model k_off
def slip_func(x,x_B,k_off_0):
    return k_off_0*np.exp((x_B*x)/(k_b*temp))

# catch-slip model k_off
def catch_slip(f,E_21,k_1rup,f_12,k_2rup,x_B):
    exp_1 = np.exp(E_21/(k_b*temp))
    exp_2 = np.exp(f/f_12)
    exp_3 = np.exp((x_B*f)/(k_b*temp))
    return (exp_1*k_1rup + exp_2*k_2rup*exp_3) / (exp_1 + exp_2)

# # slip model R^2
# def rsquared_slip(x,y,x_B,k_off_0):
#     # yfit = (x_B/(k_b*temp))*x + np.log(k_off_0)
#     yfit = k_off_0*np.exp((x_B*x)/(k_b*temp))
#     ymean=np.mean(y)
#     sstot=sum((y-ymean)**2)
#     ssres=sum((y-yfit)**2)
#     rs=1-ssres/sstot
#     return rs

# # catch-slip model R^2
# def rsquared_catch_slip(f,x,y):
#     popt,pcov = optimize.curve_fit(f,x,y,bounds=np.array([0,np.inf]))
#     residuals = y - f(x,*popt)
#     ss_res = np.sum(residuals**2)
#     ss_tot = np.sum((y-np.mean(y))**2)
#     rsq = 1 - (ss_res/ss_tot)
#     return rsq

# catch-slip fitting
# fit parameters: [E_21,k_1rup,f_12,k_2rup,x_B]
E_21_list = []
k_1rup_list = []
f_12_list = []
k_2rup_list = []
x_B_list_cs = []

for i in range(len(forces)):
    params, matrix = optimize.curve_fit(catch_slip, forces[i], 
                                                koff_avg_vals[i],
                                                bounds = np.array([0,np.inf]))
    E_21_list.append(params[0])
    k_1rup_list.append(params[1])
    f_12_list.append(params[2])
    k_2rup_list.append(params[3])
    x_B_list_cs.append(params[4])

# %% slip fitting
x_B_list_s = []
koff_0_list = []

for i in range(len(forces)):
    x_B_s, koff_0 = slip(forces[i], np.log(koff_avg_vals[i]))
    # nonspec_slip_fit = slip_func(forces, x_B_nonspec, koff_0_nonspec)
    x_B_list_s.append(x_B_s)
    koff_0_list.append(koff_0)
    
# %% distinguishing bond models
def compare_x_B(x_B_s_arr, x_B_cs_arr):
    
    # catch-slip parameters
    x_B_cs_final = []
    E_21_final = []
    k_1rup_final = []
    k_2rup_final = []
    f_12_final = []
    
    # slip model parameters
    k_off_0_final = []
    x_B_slip_final = []
    
    # colors_s = [] # colors for graphing slip
    # colors_cs = [] # colors for graphing catch-slip 
    s_peptides = [] # names of peptides (slip)
    cs_peptides = [] # names of peptides (slip)
    
    s_forces = [] # force values (slip)
    cs_forces = [] # force values (catch-slip)
    s_k_off = [] # k_off values (slip)
    cs_k_off = [] # k_off values (catch-slip)
    
    # s_exp = [] 

    if len(x_B_s_arr) != len(x_B_cs_arr):
        print('Arrays must be equal in length.')
        
    else:
        for i in range(len(x_B_cs_arr)):
            
            # if x_B is outside of a given threshold, the data is slip
            if (x_B_cs_arr[i] < 10**(-3)) or (x_B_cs_arr[i] > 10):
                x_B_slip_final.append(x_B_s_arr[i])
                k_off_0_final.append(koff_0_list[i])
                # colors_s.append(graph_colors[i])
                # s_peptides.append(peptides_list[i])
                s_forces.append(forces[i])
                s_k_off.append(koff_avg_vals[i])
                # s_exp.append(s_fit_list[i])
              
            # if x_B is inside a given threshold, the data is catch-slip
            else:
                x_B_cs_final.append(x_B_cs_arr[i])
                E_21_final.append(E_21_list[i])
                k_1rup_final.append(k_1rup_list[i])
                k_2rup_final.append(k_2rup_list[i])
                f_12_final.append(f_12_list[i])
                # colors_cs.append(graph_colors[i])
                # cs_peptides.append(peptides_list[i])
                cs_forces.append(forces[i])
                cs_k_off.append(koff_avg_vals[i])
        
        # all slip
        if (len(x_B_slip_final) != 0) and (len(x_B_cs_final) == 0):
            return k_off_0_final, x_B_slip_final, s_forces, s_k_off
    
        # all catch-slip
        elif (len(x_B_slip_final) == 0) and (len(x_B_cs_final) != 0):
            return x_B_cs_final, E_21_final, k_1rup_final, k_2rup_final, f_12_final, cs_forces, cs_k_off
     
        elif (len(x_B_slip_final) != 0) and (len(x_B_cs_final) != 0):
            return x_B_cs_final, E_21_final, k_1rup_final, k_2rup_final, f_12_final, cs_forces, cs_k_off, k_off_0_final, x_B_slip_final, s_forces, s_k_off

# results of distinguishing
lists = compare_x_B(x_B_list_s, x_B_list_cs)
if len(lists) == 4:
    koff_0_vals = lists[0]
    x_B_vals = lists[1]
    force_vals_s = lists[2]
    koff_vals_s = lists[3]
    
elif len(lists) == 7:
    x_B_vals_cs = lists[0]
    E_21_vals = lists[1]
    k_1rup_vals = lists[2]
    k_2rup_vals = lists[3]
    f_12_vals = lists[4]
    force_vals_cs = lists[5]
    koff_vals_cs = lists[6]
    
elif len(lists) == 11:
    x_B_vals_cs = lists[0]
    E_21_vals = lists[1]
    k_1rup_vals = lists[2]
    k_2rup_vals = lists[3]
    f_12_vals = lists[4]
    force_vals_cs = lists[5]
    koff_vals_cs = lists[6]
    koff_0_vals = lists[7]
    x_B_vals_s = lists[8]
    force_vals_s = lists[9]
    koff_vals_s = lists[10]
    

# %% plotting koff vs force
# slip model
if (len(lists) == 4) or (len(lists) == 11):
    plt.figure(0)
    plt.xlabel('Force (pN)')
    plt.ylabel(r'$k_{off}$')
    plt.title('Slip')
    
    for i in range(len(koff_vals_s)):
        plt.plot(force_vals_s[i], koff_vals_s[i], '.')
        
        plt.plot(force_vals_s[i],
                  slip_func(np.array(force_vals_s[i]),
                            x_B_vals_s[i],
                            koff_0_vals[i]),
                  )
        
        # avoid graphing error when only 1 trial is used
        for i in range(len(koff_error_vals)):
            for j in range(len(koff_error_vals[i])):
                if koff_error_vals[i][j] == koff_error_vals[i][j]:
                    plt.errorbar(force_vals_s[i][j],
                                 slip_func(np.array(force_vals_s[i]),
                                           x_B_vals_s[i],
                                           koff_0_vals[i]),
                                 yerr=koff_error_vals[i][j],
                                 ecolor='k',
                                 capsize=5)

if (len(lists) == 7) or (len(lists) == 11):
    plt.figure(1)
    plt.xlabel('Force (pN)')
    plt.ylabel(r'$k_{off}$')
    plt.title('Catch-slip')
    
    for i in range(len(koff_vals_cs)):
        plt.plot(force_vals_cs[i], koff_vals_cs[i], '.')
        
        plt.plot(force_vals_cs[i],
                  catch_slip(np.array(force_vals_cs[i]),
                            E_21_vals[i],
                            k_1rup_vals[i],
                            f_12_vals[i],
                            k_2rup_vals[i],
                            x_B_vals_cs[i]),
                  )
        
        # avoid graphing error when only 1 trial is used
        for i in range(len(koff_error_vals)):
            for j in range(len(koff_error_vals[i])):
                if koff_error_vals[i][j] == koff_error_vals[i][j]:
                    plt.errorbar(force_vals_cs[i][j],
                                          catch_slip(np.array(force_vals_cs[i][j]),
                                                              E_21_vals[i],
                                                              k_1rup_vals[i],
                                                              f_12_vals[i],
                                                              k_2rup_vals[i],
                                                              x_B_vals_cs[i]),
                                          yerr=koff_error_vals[i][j],
                                          ecolor='k',
                                          capsize=5)

# %% NbNT vs. k_off plot
plt.figure(2)
plt.xlabel('Force (pN)')
plt.ylabel(r'$N_b / N_T$')
for i in range(len(NbNT_vals)):
    plt.plot(forces[i], NbNT_vals[i], '.')
    plt.plot(forces[i], NbNT_vals[i])
    
for i in range(len(NbNT_error_vals)):
            for j in range(len(NbNT_error_vals[i])):
                if NbNT_error_vals[i][j] == NbNT_error_vals[i][j]:
                    plt.errorbar(forces[i][j], NbNT_vals[i][j],
                                          yerr=NbNT_error_vals[i][j],
                                          ecolor='k',
                                          capsize=5)