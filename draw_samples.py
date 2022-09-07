# -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Swiss Tropical and Public Health Institute

This code is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.

Author: Aatreyee Mimi Das

Draw samples from a Latin Hypercube for running in the model in order to calculate sensitivity indices
"""

#%%
# ----------- Import modules needed -------------

import numpy as np
from SALib.sample import saltelli
from os import chdir
import pandas as pd

#%%

# ----------- Set up problem definition for SALib module -------------

chdir('~/data/')
df = pd.read_csv('RADZEC_data.csv')
pop = df['pop']

problem = {
    'num_vars': 16,
    'names': ['betaP', 
              'betaU', 
              'betaM', 
              'inflength', 
              'tauHHP', 
              'tauHHU',
              'treatseekP', 
              'treatseekU', 
              'rdtsens',
              'followup',
              'thetaUP',
              'thetaMP',
              'thetaPU',
              'thetaMU',
              'thetaPM',
              'thetaUM'],
    'bounds': [[0.0024, 0.0072],      # point estimate: 0.0048, +/- 50%
               [0.0019, 0.0056],      # point estimate: 0.0037, +/- 50%
               [0.0027, 0.0081],    # point estimate: 0.0054, +/- 50%
               [184, 237],      # point estimate: 200 days, Sama et al, https://doi.org/10.1016/j.trstmh.2005.11.001
               [2.0, 4.8],             # point estimate: 13.1, RADZEC data, PCR positive in index HH
               [8.0, 12.7],             # point estimate: 22, RADZEC data, PCR positive in index HH
               [0.00015, 0.00044],          # point estimate: 2.9 x 10^-4, +/- 50%
               [0.00031, 0.00092],          # point estimate: 6.1 x 10^-4, +/- 50%
               [0.0, 1.0],              # point estimate: 0.34, full range
               [0.0, 1.0],              # point estimate: 0.35, full range
               [0.0016, 0.0048],           # point estimate: 0.0032, +/- 50%
               [0.0, 0.12],           # point estimate: 0.0061, Le Menach et al.
               [0.0019, 0.0058],           # point estimate: 0.0039., +/- 50%
               [0.0, 0.12],           # point estimate: 0.026, Le Menach et al.
               [0.0, 0.12*pop[0]/pop[2]],           # point estimate: 0.000057, Le Menach et al, scaled for population
               [0.0, 0.12*pop[1]/pop[2]]]           # point estimate: 0.00053, Le Menach et al, scaled for population
}
#%%
# ----------- Generate samples and save them -------------

chdir('~/sensitivity_analysis/')

n_samples2 = 2**15
param_values2 = saltelli.sample(problem, n_samples2, calc_second_order=False)
print(np.shape(param_values2))
filename = 'param_values_saltelli_' + str(n_samples2) + '.txt'
np.savetxt(filename, param_values2)