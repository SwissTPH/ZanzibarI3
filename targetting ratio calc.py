# -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Swiss Tropical and Public Health Institute
Author: Aatreyee Mimi Das

Calculating the ratio between the outputs of the targetting ratio function
presented in Chitnis et al, 2019, Malaria Journal, and the targetting ratio
observed in RADZEC data in Zanzibar

"""
#%%
# ----------- Import modules needed -------------
import numpy as np
import pandas as pd
import os
#%%
# -------- Function for targetting ratio function in Chitnis et al, 2019 ------
def targetting(a1, a2, a3, nu, pop, I):
    tau_2019 = np.exp((-1*a1*np.log(I)+(a2/nu)-(a3/nu)*np.log(I))*(pop-nu)/(pop))
    
    return tau_2019

def calc_tau(nu, scaling_factor, N, prev):
    tau_estimate = np.exp((-1*a1*np.log(prev)+(a2/nu)-(a3/nu)*np.log(prev))*(N-nu)/(N))
    tau = scaling_factor*tau_estimate
    tau_adjusted = (tau*prev*nu + 1)/(prev*(nu+1))
    return tau_adjusted

# parameters for targetting ratio function from Chitnis et al, 2019

a1 = 0.23
a2 = -1.40
a3 = 2.87

# values from RADZEC data

os.chdir('~/data/')

df = pd.read_csv('RADZEC_data.csv')
print(df)

pop = df['pop']  # Total population in each patch
I_eq = df['I_eq_neither'] #Equilibrium value of prevalence

pcr_pos_index_HH = df['pcr_pos_index_HH'] #this does not include the index case
pcr_pos_neighbour = df['pcr_pos_neighbour']
nu_HH = df['nu_HH'] - 1 # number of people investigated within index HH
nu_n = df['nu_n'] # number of people investigated in neighbouring houses
nu_tot = nu_HH+nu_n


tau_HH = pcr_pos_index_HH/I_eq  # targetting ratio within index households
tau_n_and_HH = (pcr_pos_neighbour*nu_n+pcr_pos_index_HH*nu_HH)/(I_eq*(nu_n+nu_HH))  # targetting ratio in both neighbours and index HHs, as that is what can be compared to Chitnis et al. 2019

tau_2019_HH = targetting(a1,a2,a3,nu_HH,pop,I_eq)
tau_2019_n_HH = targetting(a1,a2,a3,nu_tot,pop,I_eq)

print('Zanzibar data, HH tau: ', tau_HH, ' Zambia data, tau HH: ', tau_2019_HH)
print('Zanzibar data, HH+n tau: ', tau_n_and_HH, ' Zambia data, tau HH+n: ', tau_2019_n_HH)

print(tau_HH/tau_2019_HH)
print(tau_n_and_HH/tau_2019_n_HH) #IGNORE LAST ROW IN THE RESULTS, NO RCD ON MAINLAND TZ

scaling_factor = tau_HH/tau_2019_HH

print('Scaling factor:', scaling_factor)

#Test that the function can convert a prevalence to a targetting ratio
I_eq = 0.01
tau_calc_HH = calc_tau(nu_HH,scaling_factor,pop,I_eq)
tau_calc_n_HH = calc_tau(nu_tot,scaling_factor,pop,I_eq)
