# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 08:34:05 2020

@author: Allison
"""
import pandas as pd
import numpy as np
from scipy import stats, optimize, special
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.special as sp
import xlsxwriter


# # variables (should be inputted by user)
# a = 5  # cell radius (um)
# y = 5.025 # distance from wall to edge of cell (um)
# gamma_w = 1 # wall shear rate 
# L = 20 # bond length (nm)
# w = 800 # chamber width (um)
# b = 50 # chamber half height (um)
# mu = 6.92e-3 # (dyne/cm^2)

# %% so much user input
# empty lists to store data
track_data_files = []
spots_data_files = []

specific_site_density = []
nonspec_site_density = []

specific_force = []
nonspec_force = []
idx = 1

# coating concs
C_l_vals = []

# getting user input
# site density
check = input('Enter \"y\" if ligand site density was characterized; otherwise, enter \"n\": ')
if check.lower() == 'y':
    site_density = int(input('Enter ligand site density (sites/cm^2): '))
    specific_site_density.append(site_density)
    nonspec_site_density.append(site_density)
    
elif check.lower() == 'n':
    C_l = int(input('Enter ligand coating concentration: '))
    C_l_vals.append(C_l)
    
# flow rate
Q = int(input('Enter flow rate: ')) * (10**27 / 3600) # mL/h to pm^3/s

# Trackmate files
# nonspec_track_data = input('For flow rate = %d, enter name of \"Track statistics\" file(s) (from Trackmate) for non-specific ligand: ' % Q)
# nonspec_file_list = [val.strip() for val in nonspec_track_data.split(',')]
# nonspec_track_files = []
# for i in range(len(nonspec_file_list)):
#     if not '.csv' in nonspec_file_list[i]:
#         nonspec_file_list[i] += '.csv'
        
#     try:
#         with open(nonspec_file_list[i]) as file_open:
#             file = file_open.read()
#             nonspec_track_files.append(nonspec_file_list[i])
    
#     except FileNotFoundError:
#         print('Invalid file name.')
    
nonspec_spots_data = input('For flow rate = %d, enter name of \"Spots statistics\" file(s) (from Trackmate) for non-specific ligand: ' % Q)
nonspec_spots_list = [val.strip() for val in nonspec_spots_data.split(',')]
nonspec_spots_files = []
for i in range(len(nonspec_spots_list)):
    if not '.csv' in nonspec_spots_list[i]:
        nonspec_spots_list[i] += '.csv'
        
    try:
        with open(nonspec_spots_list[i]) as file_open:
            file = file_open.read()
            nonspec_spots_files.append(nonspec_spots_list[i])
    
    except FileNotFoundError:
        print('Invalid file name.')

# spec_track_data = input('For flow rate = %d, enter name of \"Track statistics\" file(s) (from Trackmate) for specific ligand: ' % Q)
# spec_file_list = [val.strip() for val in spec_track_data.split(',')]
# spec_track_files = []
# for i in range(len(spec_file_list)):
#     if not '.csv' in spec_file_list[i]:
#         spec_file_list[i] += '.csv'
        
#     try:
#         with open(spec_file_list[i]) as file_open:
#             file = file_open.read()
#             spec_track_files.append(spec_file_list[i])
    
#     except FileNotFoundError:
#         print('Invalid file name.')
    
spec_spots_data = input('For flow rate = %d, enter name of \"Spots statistics\" file(s) (from Trackmate) for specific ligand: ' % Q)
spec_spots_list = [val.strip() for val in spec_spots_data.split(',')]
spec_spots_files = []
for i in range(len(spec_spots_list)):
    if not '.csv' in spec_spots_list[i]:
        spec_spots_list[i] += '.csv'
        
    try:
        with open(spec_spots_list[i]) as file_open:
            file = file_open.read()
            spec_spots_files.append(spec_spots_list[i])
    
    except FileNotFoundError:
        print('Invalid file name.')
 
# %% force calculations
# assuming 1 Q value, 1 site density/coating conc
# mu = float(input('Enter viscosity (dyne-s/cm^2): ')) # check viscosity units
# a = float(input('Enter cell radius (microns): '))
# d = float(input('Enter critical distance (microns): '))
# L = float(input('Enter receptor-ligand bond length (nm): '))
# b = float(input('Enter flow chamber height (microns): '))
# b /= 2
# w = float(input('Enter flow chamber width (microns): '))
# CCD_FPS = int(input('Enter CCD FPS: '))

# %% parameters with unit conversions
mu = float(input('Enter viscosity (dyne-s/cm^2): ')) * 1e-13
a = float(input('Enter cell radius (microns): ')) * 1e6
d = float(input('Enter critical distance (microns): ')) * 1e6
L = float(input('Enter receptor-ligand bond length (nm): ')) * 1e3
b = float(input('Enter flow chamber height (microns): ')) * 10e6
b /= 2
w = float(input('Enter flow chamber width (microns): ')) * 10e6
# parameters.append([mu, a, b, L, w, d])

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

# tether force
f = Q * np.sqrt(a/(2*L)) * (1.7005 * 9*np.pi*mu*a**2 + 0.9440 * 6*np.pi*mu*a**2) / (w*b**2)
specific_force.append(f)
nonspec_force.append(f)

# %% finding non-specific lifetime
# calculate displacement of a given cell
def calc_disp(x0,x,y0,y):
        return np.sqrt((x-x0)**2+(y-y0)**2)

for i in range(len(nonspec_spots_files)):
    
    tmin_sublist = []
    spots_raw_data = pd.read_csv(nonspec_spots_files[i])
    
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
    nonspec_bond_times = []
    tc_particleID_new = []
    tc_trackID_new = []
    t_tot = 0
    
    # doing time conversion
    while i < len(tc_trackID):
        if tc_trackID[i] == tc_trackID[j]:
            if tc_frame[i]-tc_frame[k] == 1:
                t_tot += (tc_frame[i]-tc_frame[k]+6) / CCD_FPS
                if i == len(tc_trackID)-1:
                    nonspec_bond_times.append(t_tot)
                    tc_particleID_new.append(tc_particleID[k])
                    tc_trackID_new.append(tc_trackID[k])
                i += 1
                k += 1
            else:
                nonspec_bond_times.append(t_tot)
                tc_particleID_new.append(tc_particleID[k])
                tc_trackID_new.append(tc_trackID[k])
                t_tot = 0
                j = i
                i += 1
                k += 1
        else:
            nonspec_bond_times.append(t_tot)
            tc_particleID_new.append(tc_particleID[k])
            tc_trackID_new.append(tc_trackID[k])
            t_tot = 0
            j = i
            i += 1
            k += 1     
        
    t_min_avg = np.mean(nonspec_bond_times)
    tmin_sublist.append(t_min_avg)

t_min = np.mean(tmin_sublist)
print('Non-specific bond lifetime: %f' % t_min)

# %% finding specific lifetime
# track_raw_data = pd.read_csv(specific_track_data)
for i in range(len(spec_spots_files)):
    
    tmin_sublist = []
    spots_raw_data = pd.read_csv(spec_spots_files[i])
    
    # r refers to meeting criteria
    # find which particles meet stopping criteria
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
    specific_bond_times = []
    tc_particleID_new = []
    tc_trackID_new = []
    t_tot = 0
    
    # doing time conversion
    while i < len(tc_trackID):
        if tc_trackID[i] == tc_trackID[j]:
            if tc_frame[i]-tc_frame[k] == 1:
                t_tot += (tc_frame[i]-tc_frame[k]+6) / CCD_FPS
                if i == len(tc_trackID)-1:
                    specific_bond_times.append(t_tot)
                    tc_particleID_new.append(tc_particleID[k])
                    tc_trackID_new.append(tc_trackID[k])
                i += 1
                k += 1
            else:
                specific_bond_times.append(t_tot)
                tc_particleID_new.append(tc_particleID[k])
                tc_trackID_new.append(tc_trackID[k])
                t_tot = 0
                j = i
                i += 1
                k += 1
        else:
            specific_bond_times.append(t_tot)
            tc_particleID_new.append(tc_particleID[k])
            tc_trackID_new.append(tc_trackID[k])
            t_tot = 0
            j = i
            i += 1
            k += 1
            
    t_min_avg = np.mean(specific_bond_times)
    tmin_sublist.append(t_min_avg)

t_min = np.mean(specific_bond_times)
print('Specific bond lifetime: %f' % t_min)

# %% Welch's t-test on specific and non-specific bond lifetimes
# should only be comparing 1 Q, 2 ligands
t_stat, pvalue = stats.ttest_ind(nonspec_bond_times,
                                  specific_bond_times,
                                  equal_var=False)
if np.mean(specific_bond_times) <= np.mean(nonspec_bond_times):
    print('WARNING: Specific bond lifetimes were calculated to be less than non-specific bond lifetimes.')

else:
    print('Specific bond lifetimes were calculated to be greater than non-specific bond lifetimes. Your data looks reasonable!')
    
# %% writing bond lifetimes to Excel
# sheet.write(row,column,'content')
# can only write one entry at a time

# order bond lifetimes in increasing order of force
# lowest column number, lowest force value
# output force in first row
# include site density in outputted Excel file
# pandas: ordering the columns?
# Nb/NT vs koff data instead of ligand site density (which is hard to get)

# writing calculated data to excel file 


# ordering data by increasing force
# nonspec_force = sorted(nonspec_force)
# nonspec_bond_times = [time for _,time in sorted(zip(nonspec_force, nonspec_bond_times))]

# specific_force = sorted(specific_force)
# specific_bond_times = [time for _,time in sorted(zip(specific_force, specific_bond_times))]

# %% write data to Excel file
# write experimental parameters to an Excel file
variables = ['mu', 'a', 'b', 'L', 'w']

# write variable names to column 1
wb = xlsxwriter.Workbook('Module_0A_file.xlsx')
sheet0 = wb.add_worksheet('Experimental Parameters')

# headers
sheet0.write(0, 0, 'Variables')
sheet0.write(0, 1, 'Values')

for i in range(len(variables)):
    # write(row, col, content)
    sheet0.write(i+1, 0, variables[i])
    sheet0.write(i+1, 1, system_params[i])
    
# nonspec data
sheet1 = wb.add_worksheet('Nonspecific bond lifetimes')

sheet1.write(0, 0, 'Force')
for i in range(len(nonspec_force)):
    sheet1.write(0, 1, nonspec_force[i])
    
sheet1.write(1, 0, 'Ligand site density (if applicable)')
if len(nonspec_site_density) > 0:
    for i in range(len(nonspec_site_density)):
        sheet1.write(1, 1, nonspec_site_density[i])
        
sheet1.write(2, 0, 'Ligand coating conc (if applicable)')
if len(C_l_vals) > 0:
    for i in range(len(C_l_vals)):
        sheet1.write(2, 1, C_l_vals[i])
    
sheet1.write(3, 0, 'Lifetimes')
for i in range(1, len(nonspec_bond_times)):
    sheet1.write(i+2, 1, nonspec_bond_times[i])
    
# specific data
# nonspec data
sheet2 = wb.add_worksheet('Specific bond lifetimes')

sheet2.write(0, 0, 'Force')
for i in range(len(specific_force)):
    sheet2.write(0, 1, nonspec_force[i])
    
sheet2.write(1, 0, 'Ligand site density (if applicable)')
if len(nonspec_site_density) > 0:
    for i in range(len(specific_site_density)):
        sheet2.write(1, 1, specific_site_density[i])
        
sheet2.write(2, 0, 'Ligand coating conc (if applicable)')
if len(C_l_vals) > 0:
    for i in range(len(C_l_vals)):
        sheet2.write(2, 1, C_l_vals[i])
    
sheet2.write(3, 0, 'Lifetimes')
for i in range(1, len(specific_bond_times)):
    sheet2.write(i+2, 1, specific_bond_times[i])
    
wb.close()
