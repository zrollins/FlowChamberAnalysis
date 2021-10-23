# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:31:15 2020

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

# %% system params
# # variables (should be inputted by user)
# a = 5  # cell radius (um)
# y = 5.025 # distance from wall to edge of cell (um)
# gamma_w = 1 # wall shear rate 
# L = 20 # bond length (nm)
# w = 800 # chamber width (um)
# b = 50 # chamber half height (um)
# mu = 6.92e-3 # (dyne/cm^2)

# %% parameters with unit conversions
# parameters = []
mu = float(input('Enter viscosity (dyne-s/cm^2): ')) * 1e-13
a = float(input('Enter cell radius (microns): ')) * 1e6
d = float(input('Enter critical distance (microns): ')) * 1e6
L = float(input('Enter receptor-ligand bond length (nm): ')) * 1e3
b = float(input('Enter flow chamber height (microns): ')) * 1e6
b /= 2
w = float(input('Enter flow chamber width (microns): ')) * 1e6
# parameters.append([mu, a, b, L, w, d])

CCD_FPS = int(input('Enter CCD FPS: '))
m_r = int(input('Enter receptor site density (sites/um^2): '))

y = a+d

# %% loop that inputs user data
# need data with same flow rate, different site densities
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
        
        # flow rate
        Q_nc = float(input('Enter flow rate (microliter/hour): ')) 
        Q = Q_nc * (10**27 / 3600) # microliter/h to pm^3/s conversion
        
        m_l = input('For Q = %.2f, enter site densities (sites/um^2): ' % Q_nc)
        m_l_str = [val.strip() for val in m_l.split(',')]
        
        m_l_arr = np.zeros(len(m_l_str))
        for i in range(len(m_l_str)):
            m_l_arr[i] = int(m_l_str[i])
            
        site_densities.append(list(m_l_arr))
            
        # tether force
        f = Q * np.sqrt(a/(2*L)) * (1.7005*9*np.pi*mu*a**2 + 0.9440*6*np.pi*mu*a**2) / (w*b**2)
        forces.append(f)
        
        # shear stress
        tau = (3*mu*Q) / (2*w*b**2)
        shear_rate = tau / mu
        shear_rates.append(shear_rate)
        
        # Trackmate files
        track_data_sublist = []
        spots_data_sublist = []
        t_min_sublist = []
        
        for i in range(len(m_l_arr)):
            track_file_name = input('For flow rate = %.2f and site density = %d, enter name of "Track statistics" file(s) from Trackmate: ' % (Q_nc, m_l_arr[i]))
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
                
            spots_file_name = input('For flow rate = %.2f and site density = %d, enter name of "Spots in track statistics" file(s) from Trackmate: ' % (Q_nc, m_l_arr[i]))
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
            
            t_min_input = input('If applicable, for flow rate = %.2f and site density = %d, enter non-specific binding time(s). Otherwise, enter \"n\" for all trials: ' % (Q_nc, m_l_arr[i]))
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

# # for testing only!
# shear_rates = [1, 2]
# CCD_FPS = 30

for m in range(len(track_data)):
    
    koff_avg_vals_sublist = []
    koff_error_vals_sublist = []
    Nb_vals_sublist = []
    
    for n in range(len(track_data[m])):
        
        koff_avg_vals_subsub = []
        koff_error_vals_subsub = []
        Nb_vals_subsub = []
        
        # cell velocity filtering
        u_f = y*shear_rates[m]*(1-(5/16)*(a/y)**3) * 1e-6 # convert back to microns
        
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
                    
            # getting entire k_off arrays into list
            avg_k_off = np.mean(k_off_sublist)
            std_error = stats.sem(k_off_sublist)
            Nb = len(tc_trackID_unique)
            
            koff_avg_vals_subsub.append(avg_k_off)
            koff_error_vals_subsub.append(std_error)
            Nb_vals_subsub.append(Nb)
            
        koff_avg_new = np.mean(koff_avg_vals_subsub)
        Nb_avg_new = np.mean(Nb_vals_subsub)
        
        koff_avg_vals_sublist.append(koff_avg_new)
        Nb_vals_sublist.append(Nb_avg_new)
            
    koff_avg_vals.append(koff_avg_vals_sublist)
    koff_error_vals.append(koff_error_vals_sublist)
    Nb_vals.append(Nb_vals_sublist)

# %% N_T validation
def Nb_func(mrml,NT,AcKa):
    dem = AcKa*mrml
    return NT/(1+(1/dem))

# m_r = 200 # for testing purposes only!
mrml = [[3, 4, 5]]
# Nb_vals = [[83, 53], [14, 69]]

mrml_vals = []
for i in range(len(site_densities)):
    mrml = m_r * np.array(site_densities[i])
    mrml_vals.append(mrml)
    
N_T_vals = []

for i in range(len(Nb_vals)):
    params, matrix = optimize.curve_fit(Nb_func, mrml[i], Nb_vals[i],
                                        bounds=np.array([0, np.inf]))
    N_T_vals.append(params[0])
    
for i in range(len(forces)):
    print('For force = %.4f pN, N_T = %.4f.' % (forces[i], N_T_vals[i]))