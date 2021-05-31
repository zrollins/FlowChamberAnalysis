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

files = [] # only Spots in track statistics 
m_l_list = [] # 2 m_l values per Q (lists within lists)
C_l_list = []
Q_list = [] 
f_list = []

# %% system parameters from user
# assuming these values remain constant for each condition
mu = float(input('Enter viscosity (dyne/cm^2): '))
a = float(input('Enter cell radius (cm): '))
b = float(input('Enter flow chamber height (cm): '))
b /= 2
L = float(input('Enter receptor-ligand bond length (nm): '))
w = float(input('Enter flow chamber width (cm): '))
d = float(input('Enter critical distance (cm): '))
CCD_FPS = int(input('Enter CCD FPS: '))
system_params = [mu, a, b, L, w, d]

# %% input for ligand 1
while True:
    run = str(input('Enter \"y\" to input data or \"n\" to quit: '))
    
    if run.lower() == 'n':
        break
    
    elif run.lower() == 'y':
        # flow rate
        Q = float(input('Enter flow rate: '))
        Q_list.append(Q)
        
        # tether force
        f = Q * np.sqrt(a/(2*L)) * (1.7005*9*np.pi*mu*a**2 + 0.9440*6*np.pi*mu*a**2) / (w*b**2)
        f_list.append(f)

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
            
            # 1 Spots in track statistics file per site density
            files_sublist = []
            for i in range(len(m_l_sublist)):
                spots_file = input('For site density = %d, enter name of "Spots in track statistics" file from Trackmate: ' % m_l_sublist[i])
                if not '.csv' in spots_file:
                    spots_file += '.csv'
                files_sublist.append(spots_file)
        
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
                spots_file = str(input('For coating concentration = %d, enter name of "Spots in track statistics" file from Trackmate: ' % C_l_sublist[i]))
                if not '.csv' in spots_file:
                    spots_file += '.csv'
                files_sublist.append(spots_file)
            files.append(files_sublist)
                    
        else:
            print('Please enter \"c\" or \"m\".')
    
    else:
        print('Please type y/n.')
        
# %% testing that input files work
bob = pd.read_csv(files[0][1])

# %% calculating bond lifetimes
def calc_disp(x0,x,y0,y):
        return np.sqrt((x-x0)**2+(y-y0)**2)
    
bond_times = []

for q in range(len(Q_list)):
    lifetimes_per_Q = [] # lifetimes for each Q
    for n in range(len(files[q])):
        spots_raw_data = pd.read_csv(files[q][n])
        
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
                
        lifetimes_per_Q.append(lifetimes_per_ml)
    bond_times.append(lifetimes_per_Q)
    
# %% comparing lifetimes for 2 site densities/coating concs (for each Q)

lifetime_bin_vals = []

threshold = 0.05 # for p < 0.05
t_stats = []
p_values = [] 

interval = 0.3 # size of each bin (kinda like size of bins in a histogram)

for q in range(len(Q_list)):

    lifetime_bin_vals_sublist = []
    plt.figure(q)
    for f in range(len(files[q])):
        
        # bin lifetimes
        lifetime_series = pd.Series(bond_times[q][f])
        num_bins = int((max(lifetime_series) - min(lifetime_series)) / interval)
        bins = lifetime_series.value_counts(normalize=True,sort=False,bins=num_bins)
        lifetime_bin_vals_sublist.append(bins.values)
            
        # plotting the distributions
        lifetimes_plt = np.linspace(min(lifetime_series), max(lifetime_series), len(bins.values))
        plt.title('Bond lifetimes')
        plt.plot(lifetimes_plt, bins.values, '.', label='Ligand %d' % f)
        plt.legend()
        
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