# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 08:34:05 2020

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
# make sure to check that I'm in the right directory when testing
# refer to module_1.py
# bob = pd.read_csv('control_spots2.csv')

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
# parameters = []

# experimental parameters
# assume these remain constant
# mu = float(input('Enter viscosity (dyne-s/cm^2): ')) # check viscosity units
# a = float(input('Enter cell radius (microns): '))
# d = float(input('Enter critical distance (microns): '))
# L = float(input('Enter receptor-ligand bond length (nm): '))
# b = float(input('Enter flow chamber height (microns): '))
# b /= 2
# w = float(input('Enter flow chamber width (microns): '))
# parameters.append([mu, a, b, L, w, d])

# CCD_FPS = int(input('Enter CCD FPS: '))
# m_r = int(input('Enter receptor site density (sites/um^2): '))

# y = a+d

# # inter/extrapolation (Table 3 in Goldman & Brenner pt.2)
# da = np.array([0, 10e-8, 10e-7, 10e-6, 
#                10e-5, 10e-4, 10e-3, 0.003202])
# speed_constants = np.array([0.5676, 0.5556, 0.5539, 0.5518, 
#                            0.5488, 0.5445, 0.5375, 0.5315])
# fit_func = interpolate.interp1d(da, speed_constants,
#                                 fill_value='extrapolate')

# speed_const = fit_func(d/a)

# %% parameters with unit conversions
mu = float(input('Enter viscosity (dyne-s/cm^2): ')) * 1e-13
a = float(input('Enter cell radius (microns): ')) * 1e6
d = float(input('Enter critical distance (microns): ')) * 1e6
L = float(input('Enter receptor-ligand bond length (nm): ')) * 1e3
b = float(input('Enter flow chamber height (microns): ')) * 1e6
b /= 2
w = float(input('Enter flow chamber width (microns): ')) * 1e6
# parameters.append([mu, a, b, L, w, d])

y = a+d

CCD_FPS = int(input('Enter CCD FPS: '))
m_r = int(input('Enter receptor site density (sites/um^2): ')) * 1e-12

# inter/extrapolation
da = np.array([0, 10e-8, 10e-7, 10e-6, 
                10e-5, 10e-4, 10e-3, 0.003202])
speed_constants = np.array([0.5676, 0.5556, 0.5539, 0.5518, 
                            0.5488, 0.5445, 0.5375, 0.5315])
fit_func = interpolate.interp1d(da, speed_constants,
                                fill_value='extrapolate')

speed_const = fit_func(d/a)

# %% loop that inputs user data
track_data = []
spots_data = []
site_densities = []
shear_rates = []
forces = []
N_T_vals = []
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
        Q = input('For site density = %d, enter flow rates (microliter/hour), separated by commas and without new lines: ' % m_l)
        Q_str = [val.strip() for val in Q.split(',')]
        
        Q_arr = np.zeros(len(Q_str))
        for i in range(len(Q_str)):
            # Q_arr[i] = float(Q_str[i])  
            
            # convert from microliter/h to pm^3/s
            Q_arr[i] = float(Q_str[i]) * (10**27) / 3600 
        
        Q_arr_nc = np.zeros(len(Q_str))
        for i in range(len(Q_str)):
            Q_arr_nc[i] = float(Q_str[i])
            
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
        N_T_sublist = []
        
        for i in range(len(Q_arr)):
            track_file_name = input('For flow rate = %.2f and site density = %d, enter name of "Track statistics" file (s) from Trackmate: ' % (Q_arr_nc[i], m_l))
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
                
            spots_file_name = input('For flow rate = %.2f and site density = %d, enter name of "Spots in track statistics" file (s) from Trackmate: ' % (Q_arr_nc[i], m_l))
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
            
            t_min_input = input('If applicable, for flow rate = %.2f and site density = %d, enter non-specific binding time(s). Otherwise, enter \"n\" for all trials: ' % (Q_arr_nc[i], m_l))
            t_min_str = [val.strip() for val in t_min_input.split(',')]
            t_min_subsub = []
            
            for j in range(len(t_min_str)):
                if t_min_str[j] == 'n':
                    t_min_subsub.append(0.2)
                
                else:
                    t_min_float = float(t_min_str[i])
                    t_min_subsub.append(t_min_float)
                    
            t_min_sublist.append(t_min_subsub)
            
            N_T_input = input('For flow rate = %.2f and site density = %d, enter N_T values(s): ' % (Q_arr_nc[i], m_l))
            N_T_str = [val.strip() for val in N_T_input.split(',')]
            N_T_subsub = []
            
            for j in range(len(N_T_str)):
                N_T_float = float(N_T_str[j])
                N_T_subsub.append(N_T_float)
                
            N_T_sublist.append(N_T_subsub)
        
        track_data.append(track_data_sublist)
        spots_data.append(spots_data_sublist)
        t_min_vals.append(t_min_sublist)
        N_T_vals.append(N_T_sublist)

    else:
        print('Please enter \"y\" or \"n\".')

# %% calculating k_off
koff_avg_vals = []
koff_error_vals = []
Nb_vals = []
NbNT_vals = []
NbNT_error_vals = []
u_f_vals = []
U_cell_vals = []
U_hd_vals = []

for m in range(len(track_data)):
    
    koff_avg_vals_sublist = []
    koff_error_vals_sublist = []
    Nb_vals_sublist = []
    NbNT_error_vals_sublist = []
    NbNT_vals_sublist = []
    u_f_sublist = []
    
    U_hd_sublist = []
    U_cell_sublist = []
    
    u_f_actual_sublist = []
    
    for n in range(len(track_data[m])):
        
        koff_avg_vals_subsub = []
        koff_error_vals_subsub = []
        NbNT_vals_subsub = []
        
        U_hd_sublist2 = []
        U_cell_sublist2 = []
        
        # cell velocity filtering
        u_f = y*shear_rates[m][n]*(1-(5/16)*(a/y)**3) * 1e-6 # convert back to microns
        u_f_sublist.append(u_f)
        
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
            
            # i = 0
            i = 1
            i_max = len(trackID_new)
            j = 0
            
            # calculate displacement of a given cell
            def calc_disp(x0,x,y0,y):
                return np.sqrt((x-x0)**2+(y-y0)**2)
            
            # filter using stopping criteria
            tmin_frames = t_min_vals[m][n][p] * CCD_FPS
            while i < i_max-1:
                disp1 = calc_disp(x_new[i],x_new[j],y_new[i],y_new[j])
                if disp1 <= 1:
                    # i += 1
                    # disp2 = calc_disp(x_new[i],x_new[j],y_new[i],y_new[j])
                    if i-j > tmin_frames:
                        r_particleID.append(particleID_new[i])
                        r_trackID.append(trackID_new[i])
                        r_pos_x.append(x_new[i])
                        r_pos_y.append(y_new[i])
                        r_frame.append(frame_new[i])
                    i += 1
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
            # tc_particleID_new = []
            tc_trackID_new = []
            t_tot = 0
            
            # # testing
            # t_min_frames = 2
            
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
            
            # getting times when cell is not bound (k+ validation)
            durations_raw = tracks_raw_data['TRACK_DURATION'] # units: seconds
            disp_raw = tracks_raw_data['TRACK_DISPLACEMENT']
            trackID_raw = tracks_raw_data['TRACK_ID']
            track_speeds_raw = tracks_raw_data['TRACK_MEAN_SPEED']
            
            time_moving_vals = []
            durations_sublist = []
            disp_unique = []
            
            # U_cell_subsub = []
            # U_hd_subsub = []
            
            for q in range(len(tc_trackID_unique)):
                
                for i in range(len(durations_raw)):
                    if trackID_raw[i] == tc_trackID_unique[q]:
                        U_cell = track_speeds_raw[q]
                        durations_sublist.append(durations_raw[i])
                        time_moving = (durations_raw[i] - t_total_unique[q])
                        time_moving_vals.append(time_moving)
                        disp_unique.append(disp_raw[i])
                        
                        U_hd = disp_raw[i] / time_moving
                        
                # U_hd_subsub.append(U_hd)
                # U_cell_subsub.append(U_cell)
                
                # U_hd_sublist.append(U_hd)
                # U_cell_sublist.append(U_cell)
                
                U_hd_sublist2.append(U_hd)
                U_cell_sublist2.append(U_cell)
                
            # U_hd_sublist2.append(U_hd_subsub)
            # U_cell_sublist2.append(U_cell_subsub)
            
            # computing k_off = 1/lifetimes
            k_off_sublist = []
            for q in range(len(t_total_unique)):
                if t_total_unique[q] != 0:
                    k_off_sublist.append(1/t_total_unique[q])
                else:
                    k_off_sublist.append(0) 
                    # if lifetime = 0, just append 0 as a placeholder
                    # there probably shouldn't be lifetimes of 0 anyway?
                    # maybe add something in stopping criteria loop to prevent lifetime = 0
            
            # getting entire k_off arrays into list
            avg_k_off = np.mean(k_off_sublist)
            std_error = stats.sem(k_off_sublist)
            Nb = len(tc_trackID_unique)
            NT = N_T_vals[m][n][p]
            
            koff_avg_vals_subsub.append(avg_k_off)
            koff_error_vals_subsub.append(std_error)
            NbNT_vals_subsub.append(Nb/NT)
          
        # # take averages of U_hd
        # U_hd_averages = []
        # for i in range(len(U_hd_sublist2[0])):
        #     U_hd_avg_list = [sublist[i] for sublist in U_hd_sublist2]
        #     U_hd_averages.append(np.mean(U_hd_avg_list))
        
        # U_hd_sublist.append(U_hd_averages)
        
        # # take averages of U_cell
        # U_cell_averages = []
        # for i in range(len(U_cell_sublist2[0])):
        #     U_cell_avg_list = [sublist[i] for sublist in U_cell_sublist2]
        #     U_cell_averages.append(np.mean(U_cell_avg_list))
        
        # U_cell_sublist.append(U_cell_averages)
        
        koff_avg_new = np.mean(koff_avg_vals_subsub)
        koff_avg_error = stats.sem(koff_avg_vals_subsub)
        NbNT_avg_new = np.mean(NbNT_vals_subsub)
        NbNT_avg_error = stats.sem(NbNT_vals_subsub)
        
        koff_avg_vals_sublist.append(koff_avg_new)
        koff_error_vals_sublist.append(koff_avg_error)
        NbNT_error_vals_sublist.append(NbNT_avg_error)
        NbNT_vals_sublist.append(NbNT_avg_new)
        
        U_hd_sublist.append(U_hd_sublist2)
        U_cell_sublist.append(U_cell_sublist2)
        
    koff_avg_vals.append(koff_avg_vals_sublist)
    koff_error_vals.append(koff_error_vals_sublist)
    NbNT_error_vals.append(NbNT_error_vals_sublist)
    NbNT_vals.append(NbNT_vals_sublist)
    
    U_hd_vals.append(U_hd_sublist)
    U_cell_vals.append(U_cell_sublist)
    u_f_vals.append(u_f_sublist)
    
# %% k_off fitting
k_b = 0.0138
temp = 310 # Kelvin

# for testing only
# forces = [[4, 6, 9]]

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

# slip fitting
x_B_list_s = []
koff_0_list = []

for i in range(len(forces)):
    x_B_s, koff_0 = slip(forces[i], np.log(koff_avg_vals[i]))
    # nonspec_slip_fit = slip_func(forces, x_B_nonspec, koff_0_nonspec)
    x_B_list_s.append(x_B_s)
    koff_0_list.append(koff_0)
    
# distinguishing bond models
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
    # s_peptides = [] # names of peptides (slip)
    # cs_peptides = [] # names of peptides (slip)
    
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

# plotting koff vs force
# slip model
if (len(lists) == 4) or (len(lists) == 11):
    plt.figure(0)
    plt.xlabel('Force (pN)')
    plt.ylabel(r'$k_{off} (s^{-1})$')
    plt.title('Slip')
    
    for i in range(len(koff_vals_s)):
        plt.plot(force_vals_s[i], koff_vals_s[i], '.',
                 label=r'$m_l = %d$' % site_densities[i])
        
        plt.plot(force_vals_s[i],
                  slip_func(np.array(force_vals_s[i]),
                            x_B_vals_s[i],
                            koff_0_vals[i]),
                  label=r'$Best-fit curve for m_l = %d$' % site_densities[i])
        
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
    plt.legend()

if (len(lists) == 7) or (len(lists) == 11):
    plt.figure(1)
    plt.xlabel('Force (pN)')
    plt.ylabel(r'$k_{off} (s^{-1})$')
    plt.title('Catch-slip')
    
    for i in range(len(koff_vals_cs)):
        plt.plot(force_vals_cs[i], koff_vals_cs[i], '.',
                 label=r'$m_l = %d$' % site_densities[i])
        
        plt.plot(force_vals_cs[i],
                  catch_slip(np.array(force_vals_cs[i]),
                            E_21_vals[i],
                            k_1rup_vals[i],
                            f_12_vals[i],
                            k_2rup_vals[i],
                            x_B_vals_cs[i]),
                  label=r'Best-fit curve for $m_l = %d$' % site_densities[i])
        
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
                                          capsize=5, 
                                          fmt='none')
    plt.legend()
    
# %% calculating k+
def k_plus_func(NbNT, k_off, m_l):
    num = NbNT*k_off
    dem = m_l*(1-NbNT)
    return num/dem

kplus_vals = []
for i in range(len(NbNT_vals)):
    kplus = k_plus_func(np.array(NbNT_vals[i]),
                        koff_avg_vals[i],
                        site_densities[i])
    # error = stats.sem(kplus)
    kplus_vals.append(kplus)
    
# this is just for testing functionality
# be sure to delete this when we have functional data!
# kplus_vals = np.abs(kplus_vals, dtype=object)

kplus_vals_avg = []
kplus_errors = []
for i in range(len(kplus_vals[0])):
    avg_sublist = [sublist[i] for sublist in kplus_vals]
    kplus_vals_avg.append(np.mean(avg_sublist))
    kplus_errors.append(stats.sem(avg_sublist))
    
# %% calculating k_in
# # variables for Part 2 (tabulated in Cheung)
# D = 5*10**(-10) # receptor diffusivity (cm^2/s)
# a_hammer = 7.4*10**(-8) # reactive radius (cm)

# unit conversions to pm
# D = float(input('Enter diffusivity (cm^2/s): ')) * (1e10)**2 # convert to pm^2/s
# alpha = float(input('Enter reactive radius (cm): ')) * 1e10 # convert to pm

D = float(input('Enter diffusivity (cm^2/s): ')) * (1e4)**2 # convert to micron^2/s
alpha = float(input('Enter reactive radius (cm): ')) * 1e4 # convert to micron

# functions for k_in calculation
def Pe_func(u_f_arr,alpha,D):
    Pe = np.zeros(len(u_f_arr))
    for i in range(len(u_f_arr)):
        u_f_arr *= speed_const 
        Pe[i] = u_f_arr[i]*alpha/D
    
    return Pe

def Nu_func(Pe):
    I_0 = special.iv(0,Pe/2)
    K_0 = special.kv(0,Pe/2)
    summation_Nu = 0

    for n in range(1,100):
        I_n = special.iv(n,Pe/2)
        K_n = special.kv(n,Pe/2)
        summation_Nu += (-1)**n*(I_n/K_n)
    Nu = 2*(I_0/K_0 + 2*summation_Nu)
    
    return Nu

def lambda_func(Pe):
    I_0 = special.iv(0,Pe/2)
    I_1 = special.iv(1,Pe/2)
    fraction_1 = -I_1**3/I_0
    summation_lambda = 0
    
    for n in range(1,100): # max: 74, otherwise: runtime warning
        I_smol = special.iv(n-1,Pe/2)
        I_large = special.iv(n+1,Pe/2)
        I_n = special.iv(n,Pe/2)
        num = I_smol*I_large*(I_smol+I_large)
        dem = I_n
        summation_lambda += ((-1)**(n+1))*(num/dem)

    cap_lambda = (1/Pe)*(fraction_1 + summation_lambda)
    
    return cap_lambda

def k_plus_Pe(Pe,delta):
    Nu = Nu_func(Pe)
    cap_lambda = lambda_func(Pe)
    P = (cap_lambda*delta)/(1+cap_lambda*delta)
    return np.pi*D*Nu*P

def k_in_func(alpha,D,delta):
    return (delta*D)/alpha**2

Pe_vals = []
for i in range(len(u_f_vals)):
    Pe = Pe_func(u_f_vals[i], alpha, D)
    Pe_vals.append(list(Pe))
    
# %%
# for testing!
Pe_vals = [[0.72, 3]] # testing only!

Nu_vals = []
for i in range(len(Pe_vals)):
    Nu = Nu_func(np.array(Pe_vals[i]))
    Nu_vals.append(list(Nu))
    
lambda_vals = []
for i in range(len(Pe_vals)):
    lamb = lambda_func(np.array(Pe_vals[i]))
    lambda_vals.append(list(lamb))

# finding Damkohler number
delta_vals = []
for i in range(len(Pe_vals)):
    delta, matrix = optimize.curve_fit(k_plus_Pe, Pe_vals[i], kplus_vals[i])
    delta_vals.append(delta[0])
    
kin_vals = []
for i in range(len(delta_vals)):
    kin = k_in_func(alpha, D, delta_vals[i])
    kin_vals.append(kin)
    
kin_avg = np.mean(kin_vals)

# %% k+ validation
def kplus_star_func(U_hd,U_cell,k_off,m_l):
    ratio = np.array(U_hd) / np.array(U_cell)
    # ratio = U_hd / U_cell
    num = ratio*k_off*(1-(1/ratio))
    return num/m_l

kplus_star_vals = []
for i in range(len(koff_avg_vals)):
    kplus_star_sub = []
    for j in range(len(koff_avg_vals[i])):
        kplus_star = kplus_star_func(U_hd_vals[i][j], U_cell_vals[i][j], 
                                     koff_avg_vals[i][j],
                                     site_densities[i])
        kplus_star_sub.append(kplus_star)
    kplus_star_vals.append(kplus_star_sub)
    # kplus_star_vals.append(kplus_star)
    
    # U_hd_avg = np.mean(U_hd_vals[i])
    # U_cell_avg = np.mean(U_cell_vals[i])
    
kplus_star_avg_vals = []
for i in range(len(kplus_star_vals)):
    kplus_star_avg_sublist = []
    
    for j in range(len(kplus_star_vals[i])):
        kplus_star_avg = np.mean(kplus_star_vals[i][j])
        kplus_star_avg_sublist.append(kplus_star_avg)
        
    kplus_star_avg_vals.append(kplus_star_avg_sublist)

kplus_star_avg2 = []
for i in range(len(kplus_star_avg_vals[0])):
    avg_sublist = [sublist[i] for sublist in kplus_star_avg_vals]
    kplus_star_avg2.append(np.mean(avg_sublist))
    
# %% k+ fitting
# P_fit = []
# for i in range(len(lambda_vals)):
#     P_sublist = kin_avg / (kin_avg + (1/np.array(lambda_vals[i])))
#     P_fit.append(P_sublist)
    
# kon_fit = []
# for i in range(len(P_fit)):
#     kon_sublist = np.pi * D * np.array(Nu_vals[i]) * np.array(P_fit[i])
#     kon_fit.append(kon_sublist)
    
def kplus_fit_func(piDNu, A_c, P):
    return A_c * piDNu * P
        
kplus_fit_vals = []
P_fit_vals = []
for i in range(len(Nu_vals)):
    fit_params = optimize.curve_fit(kplus_fit_func, 
                                       np.pi*D*np.array(Nu_vals[i]), 
                                       kplus_vals[i])[0]
    
    # kplus_fit.append(kplus_fit_sublist)
    A_c = fit_params[0]
    
    P = fit_params[1]
    
    kplus_fit = kplus_fit_func(np.pi*D*np.array(Nu_vals[i]),
                               A_c, P)
    
    kplus_fit_vals.append(kplus_fit)
    P_fit_vals.append(P)
    
kplus_fit_avg_vals = []
for i in range(len(kplus_fit_vals[0])):
    avg_sublist = [sublist[i] for sublist in kplus_fit_vals]
    kplus_fit_avg_vals.append(np.mean(avg_sublist))
    
def kin_from_P(P, lambda_vals):
    dem = np.array(lambda_vals) * (1 - 1/P)
    return -1/dem

kin_fit_vals = []
for i in range(len(lambda_vals)):
    kin_fit = kin_from_P(P_fit_vals[i], lambda_vals[i])
    kin_fit_vals.append(np.mean(kin_fit))
    
# %% NbNT fitting
NbNT_fit = []
for i in range(len(kplus_fit_vals)):
    NbNT_fit_sublist = (np.array(kplus_fit_vals[i]) * site_densities[i]) / (np.array(kplus_fit_vals[i]) * site_densities[i] + np.array(koff_avg_vals[i]))
    NbNT_fit.append(NbNT_fit_sublist)
    
# NbNT vs. k_off plot
plt.figure(2)
plt.xlabel('Force (pN)')
plt.ylabel(r'$N_b / N_T$')
for i in range(len(NbNT_vals)):
    plt.plot(forces[i], NbNT_vals[i], '.', 
             label=r'$m_l = %d$' % site_densities[i])
    plt.plot(forces[i], NbNT_fit[i],
             label=r'Best-fit curve for $m_l = %d$' % site_densities[i])
    
for i in range(len(NbNT_error_vals)):
            for j in range(len(NbNT_error_vals[i])):
                if NbNT_error_vals[i][j] == NbNT_error_vals[i][j]:
                    plt.errorbar(forces[i][j], NbNT_vals[i][j],
                                          yerr=NbNT_error_vals[i][j],
                                          ecolor='k', fmt='none',
                                          capsize=5)
plt.legend()
    
# k+ vs force plot
plt.figure(3)
plt.xlabel('Force (pN)')
plt.ylabel(r'$k_+$')
plt.plot(forces[0], kplus_vals_avg, '.')
plt.plot(forces[0], kplus_fit_avg_vals)
plt.errorbar(forces[0], kplus_vals_avg, yerr=kplus_errors, ecolor='k',
             capsize=5, fmt='none')