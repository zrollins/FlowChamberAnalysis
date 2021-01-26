# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 08:34:05 2020

@author: Allison
"""
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.special as sp
import xlwt

# variables
a = 5  # cell radius (um)
y = 5.025 # distance from wall to edge of cell (um)
gamma_w = 1 # wall shear rate 
L = 20 # bond length (nm)
w = 800 # chamber width (um)
b = 50 # chamber half height (um)
mu = 6.92e-3 # (dyne/cm^2)

track_data = []
spots_data = []
control_site_density = []
siinfekl_site_density = []
control_force = []
siinfekl_force = []
idx = 1

while True:
    run = str(input('Enter \"y\" to input data or \"n\" to quit: '))
    
    if run.lower() == 'n':
        break
    
    elif run.lower() == 'y':
        # site density
        site_density = int(input('Enter site density (sites/um^2): '))
        control_site_density.append(site_density)
        siinfekl_site_density.append(site_density)
                        
        # force
        tether_force = int(input('Enter tether force (pN): '))
        control_force.append(tether_force)
        siinfekl_force.append(tether_force)
        
        # files
        file_name = str(input('Enter file names separated by a comma (and without pressing enter) in this order: \n1. CONTROL spots data, \n2. CONTROL tracks data, \n3. SIINFEKL spots data, \n4. SIINFEKL tracks data\n\n'))
        file_list = [x.strip() for x in file_name.split(',')]
        
        for i in range(len(file_list)):
            if not '.csv' in file_list[i]:
                file_list[i] += '.csv'
        
        try:
            for j in range(len(file_list)):
                with open(file_list[j]) as file_open:
                    file = file_open.read()
                    if j % 2 == 0:
                        spots_data.append(file_list[j])
                    else:
                        track_data.append(file_list[j])
            idx += 1
        
        except FileNotFoundError:
            print('Invalid file name.')
                 
    else:
        print('Please enter \"y\" or \"n\".')

c_tracks = track_data
c_spots = spots_data

idx = 0 # if idx is odd, control; if idx is even, siinfekl
idx2 = 1 # idx2 += 1 every other iteration since there is CONTROL & SIINFEKL video for each condition

control_bond_times = []
control_k_off = []
control_k_off_avg = []
control_error = []
control_Nb = []
control_NbNT = []

siinfekl_bond_times = []
siinfekl_k_off = []
siinfekl_k_off_avg = []
siinfekl_error = []
siinfekl_Nb = []
siinfekl_NbNT = []

for i in range(len(c_tracks)):
    # figure out how to lose dependence on order of c_tracks, c_spots
    c_tracks_data = pd.read_csv(c_tracks[i])
    c_spots_data = pd.read_csv(c_spots[i])
    
    # this part is for printing k_off lists with appropriate category
    if 'control' in c_tracks[i]:
        data_category = 'control'
    elif 'siinfekl' in c_tracks[i]:
        data_category = 'siinfekl'
    
    # cell velocity filtering
    u_f = y*gamma_w*(1-(5/16)*(a/y)**3)
    
    filtered_speeds = c_tracks_data[c_tracks_data['TRACK_MEAN_SPEED'] < np.absolute(u_f)]
    filtered_tracks_list = list(filtered_speeds['TRACK_ID'])
    
    # only collect data from control_spots_in_tracks.csv for certain velocities
    better_tracks = []
    
    # obtaining track ID's present in both spots and track stats spreadsheets
    # obtaining track ID's present in both spots and track stats spreadsheets
    trackID = c_spots_data['TRACK_ID']
    particleID = c_spots_data['ID']
    x_pos = c_spots_data['POSITION_X']
    y_pos = c_spots_data['POSITION_Y']
    frame = c_spots_data['FRAME']
    for i in range(len(filtered_tracks_list)): 
        for j in range(len(trackID)):
            if trackID[j] == filtered_tracks_list[i]:
                if j != 0:
                    if trackID[j-1] != trackID[j]:
                        better_tracks.append(trackID[j])
                else:
                    better_tracks.append(trackID[j])
                    
    # new lists to categorize info
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
    # iterate through columns in speed_filtering_final
    # find which particles meet stopping criteria
    r_pos_x = []
    r_pos_y = []
    r_trackID = []
    r_particleID = []
    r_frame = []
    
    i = 0
    i_max = len(trackID_new)
    j = 0
    
    # calculate displacement
    def calc_disp(x0,x,y0,y):
        return np.sqrt((x-x0)**2+(y-y0)**2)
    
    # the code
    while i < i_max-1:
        disp1 = calc_disp(x_new[i+1],x_new[j],y_new[i+1],y_new[j])
        if disp1 <= 1:
            i += 1
            disp2 = calc_disp(x_new[i],x_new[j],y_new[i],y_new[j])
            if i-j > 6:
                r_particleID.append(particleID_new[i])
                r_trackID.append(trackID_new[i])
                r_pos_x.append(x_new[i])
                r_pos_y.append(y_new[i])
                r_frame.append(frame_new[i])
        else:
            i += 1
            j = i-1
    
    # compile stopping criteria tracks data into pandas
    # r_data = {}
    # r_data['Particle ID'] = np.array(r_particleID)
    # r_data['Track ID'] = np.array(r_trackID)
    # r_data['X position'] = np.array(r_pos_x)
    # r_data['Y position'] = np.array(r_pos_y)
    # r_data['Frame'] = np.array(r_frame)
    # r_panda = pd.DataFrame(r_data)I
    
    # time conversion
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
                t_tot += (tc_frame[i]-tc_frame[k]+6)/30
                if i == len(tc_trackID)-1:
                    t_total.append(t_tot)
                    tc_particleID_new.append(tc_particleID[k])
                    tc_trackID_new.append(tc_trackID[k])
                i += 1
                k += 1
            else:
                t_total.append(t_tot)
                tc_particleID_new.append(tc_particleID[k])
                tc_trackID_new.append(tc_trackID[k])
                t_tot = 0
                j = i
                i += 1
                k += 1
        else:
            t_total.append(t_tot)
            tc_particleID_new.append(tc_particleID[k])
            tc_trackID_new.append(tc_trackID[k])
            t_tot = 0
            j = i
            i += 1
            k += 1
    
    # getting all bond lifetimes into list
    if idx % 2 == 0:
        control_bond_times.append(t_total)
    else:
        siinfekl_bond_times.append(t_total)
        
    # computing k_off
    k_off = []
    for n in range(len(t_total)):
        if t_total[n] != 0:
            k_off.append(1/t_total[n])
        else:
            k_off.append(0)
    
    # getting entire k_off arrays into list
    if idx % 2 == 0:
        control_k_off.append(k_off)
    else:
        siinfekl_k_off.append(k_off)
        
    avg_k_off = np.mean(k_off)
    std_error = stats.sem(k_off)
    NT = len(r_particleID)
    Nb = len(k_off)
    
    if idx % 2 == 0:
        control_k_off_avg.append(avg_k_off)
        control_error.append(std_error)
        control_Nb.append(Nb)
        control_NbNT.append(Nb/NT)
    else:
        siinfekl_k_off_avg.append(avg_k_off)
        siinfekl_error.append(std_error)
        siinfekl_Nb.append(Nb)
        siinfekl_NbNT.append(Nb/NT)
        
    if idx % 2 == 1:
        idx2 += 1
        
    idx += 1

# control
# control_data = {}
# control_data['CONTROL force'] = control_force
# control_data['Site density'] = control_site_density
# control_data['k_off'] = control_k_off
# control_data['Std error'] = control_error
# control_data['Nb/NT'] = control_NbNT
# control = pd.DataFrame(control_data)
# print(control)
# print()

# # siinfekl
# siinfekl_data = {}
# siinfekl_data['SIINFEKL force'] = siinfekl_force
# siinfekl_data['Site density'] = siinfekl_site_density
# siinfekl_data['k_off'] = siinfekl_k_off
# siinfekl_data['Std error'] = siinfekl_error
# siinfekl_data['Nb/NT'] = siinfekl_NbNT
# siinfekl = pd.DataFrame(siinfekl_data)
# print(siinfekl)
# print()

# writing bond lifetimes to Excel
# sheet.write(row,column,'content')
# can only write one entry at a time

# order bond lifetimes in increasing order of force
# lowest column number, lowest force value
# output force in first row
# include site density in outputted Excel file
# pandas: ordering the columns?
# Nb/NT vs koff data instead of ligand site density (which is hard to get)

# ordering
control_force = sorted(control_force)
control_bond_times = [time for _,time in sorted(zip(control_force,control_bond_times))]
control_k_off = [koff for _,koff in sorted(zip(control_force,control_k_off))]

siinfekl_force = sorted(siinfekl_force)
siinfekl_bond_times = [time for _,time in sorted(zip(siinfekl_force,siinfekl_bond_times))]
siinfekl_k_off = [koff for _,koff in sorted(zip(siinfekl_force,siinfekl_k_off))]

wb = xlwt.Workbook()
# CONTROL data
# bond lifetimes
sheet1 = wb.add_sheet('CONTROL bond lifetimes')

sheet1.write(0,0,'Force')
for i in range(len(control_force)):
    sheet1.write(0,i+1,control_force[i])
    
sheet1.write(1,0,'Ligand site density (if applicable)')
if len(control_site_density) > 0:
    for i in range(len(control_site_density)):
        sheet1.write(1,i+1,control_site_density[i])
        
sheet1.write(2,0,'Nb')
for i in range(len(control_Nb)):
    sheet1.write(2,i+1,control_Nb[i])
    
sheet1.write(3,0,'Lifetimes')
row = 3
col = 1
for i in range(len(control_bond_times)):
    for j in range(len(control_bond_times[i])):
        sheet1.write(row,col,control_bond_times[i][j])
        if j == len(control_bond_times[i])-1:
            row = 3
        else:
            row += 1
    col += 1

# k_off
sheet2 = wb.add_sheet('CONTROL k_off')

sheet2.write(0,0,'Force')
for i in range(len(control_force)):
    sheet2.write(0,i+1,control_force[i])
    
sheet2.write(1,0,'Ligand site density (if applicable)')
if len(control_site_density) > 0:
    for i in range(len(control_site_density)):
        sheet2.write(1,i+1,control_site_density[i])
        
sheet2.write(2,0,'Nb')
for i in range(len(control_Nb)):
    sheet2.write(2,i+1,control_Nb[i])
    
sheet2.write(3,0,'k_off')
row = 3
col = 1
for i in range(len(control_k_off)):
    for j in range(len(control_k_off[i])):
        sheet2.write(row,col,control_k_off[i][j])
        if j == len(control_k_off[i])-1:
            row = 3
        else:
            row += 1
    col += 1
    
# SIINFEKL data
# bond lifetimes
sheet3 = wb.add_sheet('SIINFEKL bond lifetimes')

sheet3.write(0,0,'Force')
for i in range(len(siinfekl_force)):
    sheet3.write(0,i+1,control_force[i])
    
sheet3.write(1,0,'Ligand site density (if applicable)')
if len(siinfekl_site_density) > 0:
    for i in range(len(siinfekl_site_density)):
        sheet3.write(1,i+1,siinfekl_site_density[i])
        
sheet3.write(2,0,'Nb')
for i in range(len(siinfekl_Nb)):
    sheet3.write(2,i+1,siinfekl_Nb[i])
    
sheet3.write(3,0,'Lifetimes')
row = 3
col = 1
for i in range(len(siinfekl_bond_times)):
    for j in range(len(siinfekl_bond_times[i])):
        sheet3.write(row,col,siinfekl_bond_times[i][j])
        if j == len(siinfekl_bond_times[i])-1:
            row = 3
        else:
            row += 1
    col += 1

# k_off
sheet4 = wb.add_sheet('SIINFEKL k_off')

sheet4.write(0,0,'Force')
for i in range(len(siinfekl_force)):
    sheet4.write(0,i+1,siinfekl_force[i])
    
sheet4.write(1,0,'Ligand site density (if applicable)')
if len(siinfekl_site_density) > 0:
    for i in range(len(siinfekl_site_density)):
        sheet4.write(1,i+1,siinfekl_site_density[i])
        
sheet4.write(2,0,'Nb')
for i in range(len(siinfekl_Nb)):
    sheet4.write(2,i+1,siinfekl_Nb[i])
    
sheet4.write(3,0,'k_off')
row = 3
col = 1
for i in range(len(siinfekl_k_off)):
    for j in range(len(siinfekl_k_off[i])):
        sheet4.write(row,col,siinfekl_k_off[i][j])
        if j == len(siinfekl_k_off[i])-1:
            row = 3
        else:
            row += 1
    col += 1

wb.save('output_data.xls')