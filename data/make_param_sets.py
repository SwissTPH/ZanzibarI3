# -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Swiss Tropical and Public Health Institute

This code is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.

Author: Aatreyee Mimi Das

Draw samples from posterior distributions from the data for the uncertainty analysis
"""

import numpy as np
import os
import pandas as pd
from scipy.stats import beta

server = False

if server:
    os.chdir('~/data/')
else:    
    os.chdir('~/data/')

df = pd.read_csv('distributions.csv')
print(df)

rows = 1000  # number of samples to draw


#beta distribution for I_eq
I_eq = df['Positive']/df['N']
N = df['N']
N = N.astype(int)
a = df['Positive'] + 1
a = a.astype(int)
b = N - a + 2
b = b.astype(int)

I_eq_array = np.zeros(shape=(rows,2))

for i in range(rows):
    I_eq_array[i]=np.random.beta(a,b)
    
    
# beta distribution for pcr_pos
a = df['pcr_pos_index_HH_t'] + 1
b = df['pcr_pos_index_HH_b']-a+2
a = a.astype(int)
b = b.astype(int)

a2 = df['pcr_pos_neighbours_t']+1
b2 = df['pcr_pos_neighbours_b']-a2+2
a2 = a2.astype(int)
b2 = b2.astype(int)

pcr_pos_HH = np.zeros(shape=(rows,2))
pcr_pos_neighbours = np.zeros(shape=(rows,2))
for i in range(rows):
    pcr_pos_HH[i]=np.random.beta(a,b)
    pcr_pos_neighbours[i]=np.random.beta(a2,b2)


# Normal distribution for the mean for nu: N(mean, standard error)
mean_nu_HH = df['nu_index_HH_mean']
SE_nu_HH = df['nu_index_HH_se']

mean_nu_n = df['nu_neighbour_mean']
SE_nu_n = df['nu_neighbour_se']

nu_index_HH = np.zeros(shape=(rows,2))
nu_neighbour = np.zeros(shape=(rows,2))
for i in range(rows):
    nu_index_HH[i]=np.random.normal(mean_nu_HH,SE_nu_HH)
    nu_neighbour[i]=np.random.normal(mean_nu_n,SE_nu_n)
    

"""
Order of columns: I_eq, pcr_pos_index_HH,pcr_pos_neighbours,nu_index_HH,nu_neighbours
"""
par_matrix = np.hstack((I_eq_array,pcr_pos_HH,pcr_pos_neighbours,nu_index_HH,nu_neighbour))

par_matrix = pd.DataFrame(par_matrix, columns=['P_I_eq','U_I_eq',
                                               'P_pcr_pos_index_HH','U_pcr_pos_index_HH',\
                                               'P_pcr_pos_neighbours','U_pcr_pos_neighbours',\
                                               'P_nu_index_HH','U_nu_index_HH',\
                                               'P_nu_neighbours','U_nu_neighbours'])

filename = '~/data/param_' + str(rows) + '.csv'
par_matrix.to_csv(filename)


# 95% CIs
print(mean_nu_HH-2*SE_nu_HH)
print(mean_nu_HH+2*SE_nu_HH)
print(mean_nu_n-2*SE_nu_n)
print(mean_nu_n+2*SE_nu_n)

print(beta.ppf(0.025, a, b)/I_eq)
print(beta.ppf(0.975, a, b)/I_eq)
print(beta.ppf(0.025, a2, b2)/I_eq)
print(beta.ppf(0.975, a2, b2)/I_eq)
