# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 08:01:18 2020

@author: Allison
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as so
import scipy.special as sp
import scipy.stats as stats

# variables for Part 1
a = 5e-4 # cell radius (cm)
y = 5.025 # distance from wall to edge of cell (um)
L = 2e-6 # bond length (cm)
w = 0.08 # chamber width (cm)
b = 5e-3 # chamber half height (cm)
mu = 6.92e-3 # (dyne/cm^2)
gamma_w = 1 # wall shear rate 
m_r = 72 # density (um^-2)
k_b = 0.0138
temp = 310

# assume user inputted site densities are organized into list like below
site_densities = [375,140,50]
site_density_arr = np.array(site_densities)

# relevant .csv files
lawrence_1b_375 = pd.read_csv("lawrence_1b_Pselectin375.csv")
lawrence_1b_140 = pd.read_csv("lawrence_1b_Pselectin140.csv")
lawrence_1b_50 = pd.read_csv("lawrence_1b_Pselectin50.csv")

# obtain shear stress, Nb/NT from Lawrence
# NT = 183 # from Cheung
lawrence_tau_375 = lawrence_1b_375['Shear stress']
lawrence_Nb_375 = lawrence_1b_375['Nb']
avg_Nb_375 = np.mean(lawrence_Nb_375)
err_375 = stats.sem(lawrence_Nb_375)

lawrence_tau_140 = lawrence_1b_140['Shear stress']
lawrence_Nb_140 = lawrence_1b_140['Nb']
avg_Nb_140 = np.mean(lawrence_Nb_140)
err_140 = stats.sem(lawrence_Nb_140)

lawrence_tau_50 = lawrence_1b_50['Shear stress']
lawrence_Nb_50 = lawrence_1b_50['Nb']
avg_Nb_50 = np.mean(lawrence_Nb_50)
err_50 = stats.sem(lawrence_Nb_50)

# NT validation
avg_Nb = [avg_Nb_375,avg_Nb_140,avg_Nb_50]
mrml = m_r*site_density_arr

# from Cheung supp
supp_data = pd.read_csv("cheung_supp_S4.csv")
Nb_actual = supp_data['Nb']
mrml_actual = supp_data['mrml']

def Nb_func(mrml,NT,AcKa):
    dem = AcKa*mrml
    return NT/(1+(1/dem))

NT_error = np.array([err_375,err_140,err_50])
parameters, matrix = so.curve_fit(Nb_func,mrml,avg_Nb,bounds=np.array([0,np.inf]))
NT = parameters[0]

# parameters2, matrix2 = so.curve_fit(Nb_func,mrml_actual,Nb_actual,sigma=NT_error,bounds=np.array([0,np.inf]))
# NT_fit2 = parameters2[0]

plt.plot(mrml,avg_Nb,'bo')
plt.plot(mrml_actual,Nb_actual,'ro')
# plt.plot(mrml_actual,Nb_func(mrml_actual,*parameters2),'k')
# plt.plot(mrml,Nb_func(mrml,*parameters),'k')

# Nb/NT from Lawrence and validated NT
lawrence_NbNT_375 = lawrence_Nb_375/NT
lawrence_NbNT_140 = lawrence_Nb_140/NT
lawrence_NbNT_50 = lawrence_Nb_50/NT

# obtain force from Lawrence (shear stress)
def T_func(tau):
    return 273*tau

lawrence_force_375 = T_func(lawrence_tau_375)
lawrence_force_140 = T_func(lawrence_tau_140)
lawrence_force_50 = T_func(lawrence_tau_50)

# using tabulated fit parameters (Cheung) to find k_off
# cheung_constants: E_21,k_1rup,f_12,k_2rup,x_B(nm)
cheung_constants = np.array([4.456*k_b*temp,2.85,3.56,1.17,0.033])

def catch_slip(f,E_21,k_1rup,f_12,k_2rup,x_B):
    exp_1 = np.exp(E_21/(k_b*temp))
    exp_2 = np.exp(f/f_12)
    exp_3 = np.exp((x_B*f)/(k_b*temp))
    return (exp_1*k_1rup + exp_2*k_2rup*exp_3) / (exp_1 + exp_2)

lawrence_k_off_375 = catch_slip(lawrence_force_375,*cheung_constants)
lawrence_k_off_140 = catch_slip(lawrence_force_140,*cheung_constants)
lawrence_k_off_50 = catch_slip(lawrence_force_50,*cheung_constants)

# avg k_off for A_c*k_on/k_off plot
avg_k_off = np.zeros(len(lawrence_k_off_375))

for i in range(len(lawrence_k_off_375)):
    mean = np.mean([lawrence_k_off_375[i],lawrence_k_off_140[i],lawrence_k_off_50[i]])
    avg_k_off[i] = mean

# finding k_+
def k_plus_func(NbNT,k_off,m_l):
    num = NbNT*k_off
    dem = m_l*(1-NbNT)
    return num/dem

lawrence_k_plus_375 = k_plus_func(lawrence_NbNT_375,lawrence_k_off_375,m_l=375)
lawrence_k_plus_140 = k_plus_func(lawrence_NbNT_140,lawrence_k_off_140,m_l=140)
lawrence_k_plus_50 = k_plus_func(lawrence_NbNT_50,lawrence_k_off_375,m_l=50)

# finding average k_+ values and std error (including m_l=375)
avg_arr = np.zeros(len(lawrence_k_plus_375))
errors_arr = np.zeros(len(lawrence_k_plus_375))

for i in range(len(lawrence_k_plus_375)):
    mean = np.mean([lawrence_k_plus_375[i],lawrence_k_plus_140[i],lawrence_k_plus_50[i]])
    k_plus_i = np.array([lawrence_k_plus_375[i],lawrence_k_plus_140[i],lawrence_k_plus_50[i]])
    error = stats.sem(k_plus_i)
    errors_arr[i] = error
    avg_arr[i] = mean
    
k_plus_dictionary = {}
k_plus_dictionary['k_+'] = avg_arr
k_plus_dictionary['std error'] = errors_arr
k_plus_vals = pd.DataFrame(k_plus_dictionary)
print(k_plus_vals)
print()

# Ack_on stuff
Ac_kon = avg_arr/m_r
Ac_kon_koff = np.zeros(len(avg_k_off))
for i in range(len(avg_k_off)):
    Ac_kon_koff[i] = Ac_kon[i]/avg_k_off[i]
    
# k_+ cross validation
u_hd = y*gamma_w*(1-(5/16)*(a/y)**3)

    
# finding average k_+ values and std error (excluding m_l = 375)
avg_arr2 = np.zeros(len(lawrence_k_plus_140))
errors_arr2 = np.zeros(len(lawrence_k_plus_140))

for i in range(len(lawrence_k_plus_140)):
    mean = np.mean([lawrence_k_plus_140[i],lawrence_k_plus_50[i]])
    k_plus_i = np.array([lawrence_k_plus_140[i],lawrence_k_plus_50[i]])
    error = stats.sem(k_plus_i)
    errors_arr2[i] = error
    avg_arr2[i] = mean

# find average of tau values (they should be all be about equal)
tau_arr = np.zeros(len(lawrence_k_plus_375))
for i in range(len(lawrence_tau_375)):
    tau_mean = np.mean([lawrence_tau_375[i],lawrence_tau_140[i],lawrence_tau_50[i]])
    tau_arr[i] = tau_mean

# variables for Part 2 (tabulated in Cheung)
D = 0.05 # receptor diffusivity (um^2/s)
a_hammer = 7.4e-4 # reactive radius (um)

# cell velocity filtering
a_v = 5e-4 # cell radius (cm)
y_v = 5.025e-4 # distance from wall to edge of cell (cm)

# functions for Part 2n
def u_f_func(a,y,tau,mu):
    gamma_w = tau/mu
    return y*gamma_w*(1-(5/16)*(a/y)**3)*10**4 # units: um/s

def Pe_func(v_arr,a,D):
    Pe = np.zeros(len(v_arr))
    for i in range(len(v_arr)):
        Pe[i] = v_arr[i]*a/D
    return Pe

def Nu_func(Pe):
    I_0 = sp.iv(0,Pe/2)
    K_0 = sp.kv(0,Pe/2)
    summation_Nu = 0

    for n in range(1,1000):
        I_n = sp.iv(n,Pe/2)
        K_n = sp.kv(n,Pe/2)
        summation_Nu += (-1)**n*(I_n/K_n)
    Nu = 2*(I_0/K_0 + 2*summation_Nu)
    
    return Nu

def lambda_func(Pe):
    I_0 = sp.iv(0,Pe/2)
    I_1 = sp.iv(1,Pe/2)
    fraction_1 = -I_1**3/I_0
    summation_lambda = 0
    
    for n in range(1,40): # max: 74, otherwise: runtime warning
        I_smol = sp.iv(n-1,Pe/2)
        I_large = sp.iv(n+1,Pe/2)
        I_n = sp.iv(n,Pe/2)
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

def k_in_func(a,D,delta):
    return (delta*D)/a**2

u_f_375 = u_f_func(a_v,y_v,lawrence_tau_375,mu)
u_f_140 = u_f_func(a_v,y_v,lawrence_tau_140,mu)
u_f_50 = u_f_func(a_v,y_v,lawrence_tau_50,mu)

Pe_375 = Pe_func(u_f_375,a_hammer,D)
Pe_140 = Pe_func(u_f_140,a_hammer,D)
Pe_50 = Pe_func(u_f_50,a_hammer,D)

Nu_375 = Nu_func(Pe_375)
Nu_140 = Nu_func(Pe_140)
Nu_50 = Nu_func(Pe_50)

lambda_375 = lambda_func(Pe_375)
lambda_140 = lambda_func(Pe_140)
lambda_50 = lambda_func(Pe_50)

# finding delta by fitting k_+ vs. Pe
delta_375, matrix_375 = so.curve_fit(k_plus_Pe, Pe_375, lawrence_k_plus_375)
delta_140, matrix_140 = so.curve_fit(k_plus_Pe, Pe_140, lawrence_k_plus_140)
delta_50, matrix_50 = so.curve_fit(k_plus_Pe, Pe_50, lawrence_k_plus_50)
avg_delta = np.mean([delta_375,delta_140,delta_50])
delta_std_error = stats.sem([delta_375,delta_140,delta_50])

# finding k_in using delta
k_in_375 = k_in_func(a_hammer,D,delta_375)
k_in_140 = k_in_func(a_hammer,D,delta_140)
k_in_50 = k_in_func(a_hammer,D,delta_50)
avg_k_in = np.mean([k_in_375,k_in_140,k_in_50])
k_in_std_error = stats.sem([k_in_375,k_in_140,k_in_50])

# assume user inputted site densities are organized into list like below
site_densities = [375,140,50]

k_in_dictionary = {}
k_in_dictionary['site density'] = site_densities
k_in_dictionary['delta'] = [delta_375[0],delta_140[0],delta_50[0]]
k_in_dictionary['k_in'] = [k_in_375[0],k_in_140[0],k_in_50[0]]
k_in_vals = pd.DataFrame(k_in_dictionary)
print(k_in_vals)
print()
print('Mean delta:', avg_delta)
print('Delta std error:', delta_std_error)
print()
print('Mean k_in:', avg_k_in)
print('k_in std error:', k_in_std_error)

# plots
plt.figure(figsize=(10,10))

plt.subplot(4,1,1)
plt.plot(lawrence_tau_375,lawrence_NbNT_375,'bo',label='m_l = 375')
plt.plot(lawrence_tau_140,lawrence_NbNT_140,'ro',label='m_l = 140')
plt.plot(lawrence_tau_50,lawrence_NbNT_50,'go',label='m_l = 50')
plt.ylabel('Nb/NT')
plt.legend()

plt.subplot(4,1,2)
plt.plot(tau_arr,avg_arr2,'bo',label='including m_l = 375')
plt.errorbar(tau_arr,avg_arr,yerr=errors_arr,fmt='none')
plt.ylabel('k_+')
plt.legend()

plt.subplot(4,1,3)
plt.plot(tau_arr,avg_arr/(m_r*avg_k_off),'bo')
plt.ylabel('k_+/(m_r*k_off')

plt.subplot(4,1,4)
plt.plot(tau_arr,Ac_kon_koff,'bo')
plt.ylabel('A_c*k_on/k_off')
plt.xlabel('Shear stress')

plt.show()