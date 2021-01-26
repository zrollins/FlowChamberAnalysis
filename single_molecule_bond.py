# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 08:34:05 2020

@author: Allison
"""
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from xlrd import open_workbook, XLRDError

files = [] # 1 file per m_l, 2 files per Q (lists within lists)
m_l_list = [] # 2 m_l values per Q (lists within lists)
Q_list = [] 

while True:
    run = str(input('Enter \"y\" to input data or \"n\" to quit: '))
    
    if run.lower() == 'n':
        break
    
    elif run.lower() == 'y':
        # flow rate
        Q = float(input('Enter flow rate: '))
        Q_list.append(Q)

        # # site densities
        # m_l_sublist = []
        # m_l = str(input('For flow rate = %.2f, enter site densities (sites/um^2) separated by commas and without new lines: \n\n' % Q))
        # m_l_str = [val.strip() for val in m_l.split(',')]
        
        # for i in range(len(m_l_str)):
        #     m_l_sublist.append(m_l_str[i])
        # m_l_list.append(m_l_sublist)
        
        # files
        file_names = str(input('Enter file: '))
        if not '.xls' in file_names:
            file_names += '.xls'
            
        try:
            file = open_workbook(file_names)
            files.append(file_names)
        except XLRDError:
            print('Invalid file.')
            
        # file_list = [x.strip() for x in file_names.split(',')]
        
        # for i in range(len(file_list)):
        #     if not '.xls' in file_list[i]:
        #         file_list[i] += '.xls'
          
        # try:
        #     for i in range(len(file_list)):
        #         file = open_workbook(file_list[i])
        #     files.append(file_list)
        # except XLRDError:
        #     print('Invalid file.')
    
    else:
        print('Please type y/n.')
    
# validate that the forces are the same for 2 columns
# plot distributions
# rows for force/site density?
# subsequent rows for lifetimes?

lifetime_bin_vals = []
koff_bin_vals = []

# Welch
threshold = 0.05
t_stats_lifetimes = [] 
t_stats_koff = []
p_values_lifetimes = [] 
p_values_koff = []

# iterate thru lists of files at each Q value
for f in range(len(files)):
    plt.figure(figsize=(7,15))
    bins_list = [] # should only have 2 sublists, 1 for each m_l value
    
    # lifetimes
    control_bins1 = []
    siinfekl_bins1 = []
    
    # k_off
    control_bins2 = []
    siinfekl_bins2 = []
    
    # CONTROL, SIINFEKL data to compare will always be odd/even integer sheets
    # lifetimes are sheets 1&3, k_off are sheets 2&4, etc.
    
    col = 1 # iterate thru columns
    graph_idx = 1
    while True:
        # check if there is still data in each column
        placeholder = pd.read_excel(files[f],sheet_name=f,usecols=[col],skiprows=2)
        if len(placeholder) == 0:
            break
        else:
            # site densities
            m_l_empty = []
            m_l_control = pd.read_excel(files[f],sheet_name=f,usecols=[col],skiprows=1)
            for i in range(len(m_l_control)):
                m_l_empty.append(m_l_control.values[i][0])
            m_l_list.append(m_l_empty)
            
            m_l_empty = []
            m_l_siinfekl = pd.read_excel(files[f],sheet_name=f+1,usecols=[col],skiprows=1)
            for i in range(len(m_l_siinfekl)):
                m_l_empty.append(m_l_siinfekl.values[i][0])
            m_l_list.append(m_l_siinfekl)
            
            # lifetimes
            control_times = pd.read_excel(files[f],sheet_name=f,usecols=[col],skiprows=2)
            siinfekl_times = pd.read_excel(files[f],sheet_name=f+2,usecols=[col],skiprows=2)
            
            # put these values into a list to convert into series
            control_list = []
            for i in range(len(control_times)):
                control_list.append(control_times.values[i][0])
            
            siinfekl_list = []
            for i in range(len(siinfekl_times)):
                siinfekl_list.append(siinfekl_times.values[i][0])
                
            # bin values
            control_series = pd.Series(control_list)
            siinfekl_series = pd.Series(siinfekl_list)
            interval = 0.3
            num_bins1 = int(np.ceil((max(control_list)-min(control_list))/interval))
            num_bins2 = int(np.ceil((max(siinfekl_list)-min(siinfekl_list))/interval))
            bins1 = control_series.value_counts(normalize=True,sort=False,bins=num_bins1)
            control_bins1.append(bins1.values)
            bins2 = siinfekl_series.value_counts(normalize=True,sort=False,bins=num_bins2)
            siinfekl_bins1.append(bins2.values)
            
            # plotting the distributions
            control_x = np.linspace(min(control_list),max(control_list),len(bins1.values))
            siinfekl_x = np.linspace(min(siinfekl_list),max(siinfekl_list),len(bins2.values))
            plt.subplot(4*len(files),1,graph_idx)
            graph_idx += 1
            plt.title('Bond lifetimes')
            plt.plot(control_x, bins1.values, 'b', label='CONTROL')
            plt.plot(siinfekl_x, bins2.values, 'r', label='SIINFEKL')
            plt.legend()
            
            # perform Welch's t-test
            t_stat, pvalue = stats.ttest_ind(bins1.values,bins2.values,equal_var=False)
            t_stats_lifetimes.append(t_stat)
            p_values_lifetimes.append(pvalue)
                
            # if pvalue > threshold:
            #     print('Your p-value is greater than %f. Please consider retaking your data.' % threshold)
            
            # k_off
            control_koff = pd.read_excel(files[f],sheet_name=f+1,usecols=[col],skiprows=2)
            siinfekl_koff = pd.read_excel(files[f],sheet_name=f+3,usecols=[col],skiprows=2)
            
            # put these values into a list to convert into series
            control_list = []
            for i in range(len(control_koff)):
                control_list.append(control_koff.values[i][0])
            
            siinfekl_list = []
            for i in range(len(siinfekl_koff)):
                siinfekl_list.append(siinfekl_koff.values[i][0])
                
            # bin values
            control_series = pd.Series(control_list)
            siinfekl_series = pd.Series(siinfekl_list)
            interval = 0.03 # DIFFERENT THAN LIFETIMES
            num_bins1 = int(np.ceil((max(control_list)-min(control_list))/interval))
            num_bins2 = int(np.ceil((max(siinfekl_list)-min(siinfekl_list))/interval))
            bins1 = control_series.value_counts(normalize=True,sort=False,bins=num_bins1)
            control_bins1.append(bins1.values)
            bins2 = siinfekl_series.value_counts(normalize=True,sort=False,bins=num_bins2)
            siinfekl_bins1.append(bins2.values)
            
            # plotting the distributions
            control_x = np.linspace(min(control_list),max(control_list),len(bins1.values))
            siinfekl_x = np.linspace(min(siinfekl_list),max(siinfekl_list),len(bins2.values))
            plt.subplot(4*len(files),1,graph_idx)
            graph_idx += 1
            plt.title('k_off')
            plt.plot(control_x, bins1.values, 'b', label='CONTROL')
            plt.plot(siinfekl_x, bins2.values, 'r', label='SIINFEKL')
            plt.legend()
            
            # perform Welch's t-test
            t_stat, pvalue = stats.ttest_ind(bins1.values,bins2.values,equal_var=False)
            t_stats_koff.append(t_stat)
            p_values_koff.append(pvalue)
                
            # if pvalue > threshold:
            #     print('Your p-value is greater than %f. Please consider retaking your data.' % threshold)
            
        # lifetime_bin_vals.append(control_bins1, siinfekl_bins1)
        # koff_bin_vals.append(control_bins2, siinfekl_bins2)
        col += 1

plt.legend()
plt.show()
    # # iterate thru each file in a sub-list (each at 1 m_l value)
    # for i in range(len(files[f])):
    #     data = pd.read_excel(files[f][i])
    #     lifetimes = data['force %d' % column]
        
    #     # bin lifetimes at each of 2 m_l values
    #     series = pd.Series(lifetimes)
    #     interval = 0.3 # arbitrary
    #     num_bins = int(np.ceil((max(lifetimes)-min(lifetimes))/interval))
    #     bins = series.value_counts(normalize=True,sort=False,bins=num_bins)
    #     bins_list.append(bins.values)
        
    # bin_vals_list.append(bins_list)
    
    # # perform Welch's t-test on 1 Q, 2 m_l   
    # for j in range(len(bins_list)):
    #     t_stat, pvalue = stats.ttest_ind(bins_list[0],bins_list[1],equal_var=False)
    #     t_stats.append(t_stat)
    #     p_values.append(pvalue)
        
    #     if pvalue > threshold:
    #         print("Your data sucks")
            
    #     # # plt.plot(x_new,bins1.values,'b')
    #     # # plt.plot(x_new2,bins2.values,'r')
    #     # # plt.show()
                
    # column += 1
    

# # testing with just 1 file
# control_bond_times = pd.read_excel("CONTROL_bond_lifetimes.xls")
# x = control_bond_times['bond lifetimes']

# series = pd.Series(x)
# interval = 0.3 # arbitrary
# num_bins = int(np.ceil((max(x)-min(x))/interval))
# bob = series.value_counts(normalize=True,sort=False,bins=num_bins)

# create random list of numbers from 5-10 (distribution) 
# plot that using pd.cut?
# actually, use red curve from Liu et. al
# mean1 = 0.8
# mean2 = 0.4
# sd1 = 0.1
# sd2 = 0.1
# size = 100
# interval = 0.03

# red1 = np.random.normal(mean1,sd1,size)
# red2 = np.random.normal(mean2,sd2,size)

# x = pd.Series(red1)
# x2 = pd.Series(red2)

# num_bins1 = int(np.ceil((max(x)-min(x))/interval))
# num_bins2 = int(np.ceil((max(x2)-min(x2))/interval))

# bins1 = x.value_counts(normalize=True,sort=False,bins=num_bins1)
# bin_vals1 = bins1.values
# bins2 = x2.value_counts(normalize=True,sort=False,bins=num_bins)
# bin_vals2 = bins2.values

# x_new = np.linspace(0,max(x),len(bin_vals1))
# x_new2 = np.linspace(0,max(x2),len(bin_vals2))

# # Welch's t-test
# t_stat, pvalue = stats.ttest_ind(red1,red2,equal_var=False)
# threshold = 0.05
# if t_stat > threshold:
#     print("Your data sucks")

# # plt.plot(x_new,bins1.values,'b')
# # plt.plot(x_new2,bins2.values,'r')
# # plt.show()