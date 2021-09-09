# -*- coding: utf-8 -*-
"""
Created on Mon May 31 11:54:40 2021

@author: Allison
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from xlrd import open_workbook, XLRDError

# %% user input
files = [] # only Spots in track statistics 
m_l_list = [] # 2 m_l values per Q (lists within lists)
C_l_list = []
Q_list = [] 

while True:
    run = str(input('Enter \"y\" to input data or \"n\" to quit: '))
    
    if run.lower() == 'n':
        break
    
    elif run.lower() == 'y':
        # flow rate
        Q = float(input('Enter flow rate: '))
        Q_list.append(Q)

        # site densities or coating conc?
        which_one = input('Enter \"c\" to input coating concentrations or \"m\" to enter site densities: ')
        
        # characterized site densities
        if which_one == 'm':
            m_l_sublist = []
            m_l = input('For flow rate = %.2f, enter site densities (sites/um^2) separated by commas and without new lines: ' % Q)
            m_l_str = [val.strip() for val in m_l.split(',')]
            
            for i in range(len(m_l_str)):
                m_l_sublist.append(int(m_l_str[i]))
                
            m_l_list.append(m_l_sublist)
            
            # spots in track statistics file per site density
            files_sublist = []
            for i in range(len(m_l_sublist)):
                spots_file = input('For flow rate = %.2f and site density = %d, enter name of "Spots in track statistics" file(s) from Trackmate: ' % (Q, m_l_sublist[i]))
                spots_file_list = [val.strip() for val in spots_file.split(',')]
                spots_file_subsub = []
                for i in range(len(spots_file_list)):
                    if not '.csv' in spots_file_list[i]:
                        spots_file += '.csv'
                        
                    try:
                        with open(spots_file_list[i]) as file_open:
                            file = file_open.read()
                            spots_file_subsub.append(spots_file_list[i])
                    
                    except FileNotFoundError:
                        print('Invalid file name.')
                        
                    # spots_file_subsub.append(spots_file_list[i])
                
                # if not '.csv' in spots_file:
                #     spots_file += '.csv'
                
                files_sublist.append(spots_file_subsub)
                
            # files_sublist.append(spots_file_subsub)
            
            files.append(files_sublist)
        
        elif which_one == 'c':
            # coating conc
            C_l_sublist = []
            C_l = input('For flow rate = %.2f, enter coating concentrations (sites/um^2) separated by commas and without new lines: ' % Q)
            C_l_str = [val.strip() for val in C_l.split(',')]
            
            for i in range(len(C_l_str)):
                C_l_sublist.append(int(C_l_str[i]))
            C_l_list.append(C_l_sublist)
        
            # 1 Spots in track statistics file for each C_l
            files_sublist = []
            for i in range(len(C_l_sublist)):
                spots_file = input('For flow rate = %.2f and coating concentration = %d, enter name of "Spots in track statistics" file(s) from Trackmate: ' % (Q, C_l_sublist[i]))
                spots_file_list = [val.strip() for val in spots_file.split(',')]
                spots_file_subsub = []
                for i in range(len(spots_file_list)):
                    if not '.csv' in spots_file_list[i]:
                        spots_file += '.csv'
                        
                    try:
                        with open(spots_file_list[i]) as file_open:
                            file = file_open.read()
                            spots_file_subsub.append(spots_file_list[i])
                    
                    except FileNotFoundError:
                        print('Invalid file name.')
                        
                    # spots_file_subsub.append(spots_file_list[i])
                
                # if not '.csv' in spots_file:
                #     spots_file += '.csv'
                
                files_sublist.append(spots_file_subsub)
                
            # files_sublist.append(spots_file_subsub)
            
            files.append(files_sublist)
                    
        else:
            print('Please enter \"c\" or \"m\".')
    
    else:
        print('Please type y/n.')
        
# %% calculating tether force
# assume system parameters remain constant for each condition
mu = float(input('Enter viscosity (dyne-s/cm^2): ')) 
a = float(input('Enter cell radius (microns): '))
d = float(input('Enter critical distance (microns): '))
L = float(input('Enter receptor-ligand bond length (nm): '))
b = float(input('Enter flow chamber height (microns): '))
b /= 2
w = float(input('Enter flow chamber width (microns): '))
CCD_FPS = int(input('Enter CCD FPS: '))
# system_params = [mu, a, b, L, w, d]

f_list = np.array(Q_list) * np.sqrt(a/(2*L)) * (1.7005*9*np.pi*mu*a**2 + 0.9440*6*np.pi*mu*a**2) / (w*b**2)

# %% calculating bond lifetimes
def calc_disp(x0,x,y0,y):
        return np.sqrt((x-x0)**2+(y-y0)**2)
    
bond_times = []

for q in range(len(Q_list)):
    
    lifetimes_per_Q = [] # lifetimes for each Q
    
    for n in range(len(files[q])):
        
        lifetimes_subsub = []
        
        for r in range(len(files[q][n])):
            spots_raw_data = pd.read_csv(files[q][n][r])
            
            # r refers to meeting criteria
            r_pos_x = []
            r_pos_y = []
            r_trackID = []
            r_particleID = []
            r_frame = []
            
            i = 0
            j = 0
            
            x_pos = spots_raw_data['POSITION_X']
            y_pos = spots_raw_data['POSITION_Y']
            particle_ID = spots_raw_data['ID']
            track_ID = spots_raw_data['TRACK_ID']
            frames = spots_raw_data['FRAME']
            
            # number of cells in the chamber
            i_max = len(frames) 
            
            # filter using stopping criteria
            while i < i_max-1:
                disp1 = calc_disp(x_pos[i+1], x_pos[j], y_pos[i+1], y_pos[j])
                if disp1 <= 1:
                    i += 1
                    disp2 = calc_disp(x_pos[i], x_pos[j], y_pos[i], y_pos[j])
                    if i-j > 6:
                        r_particleID.append(particle_ID[i])
                        r_trackID.append(track_ID[i])
                        r_pos_x.append(x_pos[i])
                        r_pos_y.append(y_pos[i])
                        r_frame.append(frames[i])
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
            lifetimes_per_ml = []
            tc_particleID_new = []
            tc_trackID_new = []
            t_tot = 0
            
            # doing time conversion
            while i < len(tc_trackID):
                if tc_trackID[i] == tc_trackID[j]:
                    if tc_frame[i]-tc_frame[k] == 1:
                        t_tot += (tc_frame[i]-tc_frame[k]+6) / CCD_FPS
                        if i == len(tc_trackID)-1:
                            lifetimes_per_ml.append(t_tot)
                            tc_particleID_new.append(tc_particleID[k])
                            tc_trackID_new.append(tc_trackID[k])
                        i += 1
                        k += 1
                    else:
                        lifetimes_per_ml.append(t_tot)
                        tc_particleID_new.append(tc_particleID[k])
                        tc_trackID_new.append(tc_trackID[k])
                        t_tot = 0
                        j = i
                        i += 1
                        k += 1
                else:
                    lifetimes_per_ml.append(t_tot)
                    tc_particleID_new.append(tc_particleID[k])
                    tc_trackID_new.append(tc_trackID[k])
                    t_tot = 0
                    j = i
                    i += 1
                    k += 1  
            
            lifetimes_subsub.append(lifetimes_per_ml)
            
        # # take averages of U_hd
        # U_hd_averages = []
        # for i in range(len(U_hd_sublist2[0])):
        #     U_hd_avg_list = [sublist[i] for sublist in U_hd_sublist2]
        #     U_hd_averages.append(np.mean(U_hd_avg_list))
        
        # U_hd_sublist.append(U_hd_averages)
            
        avg_lifetimes = []
        for i in range(len(lifetimes_subsub[0])):
            avg_lifetimes_vals = [sublist[i] for sublist in lifetimes_subsub]
            avg_lifetimes.append(np.mean(avg_lifetimes_vals))
            
        lifetimes_per_Q.append(avg_lifetimes)
        
    bond_times.append(lifetimes_per_Q)
    
# %% comparing lifetimes for 2 site densities/coating concs (for each Q)
lifetime_bin_vals = []  

threshold = 0.05 # for p < 0.05
t_stats = []
p_values = [] 

interval = 0.3 # size of each bin (kinda like size of bins in a histogram)

# formatting plot
shapes = ['o', '+', 'o', '+'] # different site densities
colors = ['b', 'r'] # different flow rates

plt.figure(0)
plt.xlabel('Bond lifetimes (s)')
plt.ylabel('Fraction of stopping events')

for q in range(len(Q_list)):

    lifetime_bin_vals_sublist = []
    
    for f in range(len(files[q])):
        # bin lifetimes (first m_l or C_l)
        lifetime_series = pd.Series(bond_times[q][0])
        num_bins = int((max(lifetime_series) - min(lifetime_series)) / interval)
        bins = lifetime_series.value_counts(normalize=True,sort=False,bins=num_bins)
        lifetime_bin_vals_sublist.append(bins.values)
            
        # plotting the distribution
        lifetimes_plt = np.linspace(min(lifetime_series), max(lifetime_series), len(bins.values))
        
        # avoiding repeating plot labels
        if f == 0:
            # coating concentrations
            if len(C_l_list) != 0:
                plt.plot(lifetimes_plt, bins.values, '%so' % colors[q], 
                         label='$Q = %.1f, C_l = %d$' % (Q_list[q], C_l_list[q][0]))
            else:
                plt.plot(lifetimes_plt, bins.values, '%so' % colors[q], 
                         label='$Q = %.1f, m_l = %d$' % (Q_list[q], m_l_list[q][0]))
                
        else:
            plt.plot(lifetimes_plt, bins.values, '%so' % colors[q])
        
        # bin lifetimes (second m_l or C_l)
        lifetime_series = pd.Series(bond_times[q][1])
        num_bins = int((max(lifetime_series) - min(lifetime_series)) / interval)
        bins = lifetime_series.value_counts(normalize=True,sort=False,bins=num_bins)
        lifetime_bin_vals_sublist.append(bins.values)
            
        # plotting the distribution
        lifetimes_plt = np.linspace(min(lifetime_series), max(lifetime_series), len(bins.values))
        
        # avoiding repeating plot labels
        if f == 0:
            # coating concentrations
            if len(C_l_list) != 0:
                plt.plot(lifetimes_plt, bins.values, '%s+' % colors[q], 
                         label=r'$Q = %.1f, C_l = %d$' % (Q_list[q], C_l_list[q][1]))
            else:
                plt.plot(lifetimes_plt, bins.values, '%s+' % colors[q], 
                         label=r'$Q = %.1f, m_l = %d$' % (Q_list[q], m_l_list[q][1]))
                
        else:
            plt.plot(lifetimes_plt, bins.values, '%s+' % colors[q])
        
    # perform Welch's t-test
    # assume there will only be 2 site densities/coating conc per Q
    t_stat, pvalue = stats.ttest_ind(lifetime_bin_vals_sublist[0],
                                     lifetime_bin_vals_sublist[1],
                                     equal_var=False)
    
    t_stats.append(t_stat)
    p_values.append(pvalue)
        
    if pvalue < threshold:
        print('Your p-value is less than %f. Please consider retaking your data.' % threshold)
    
plt.legend()
plt.show()