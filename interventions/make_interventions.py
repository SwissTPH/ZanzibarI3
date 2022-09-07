# -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Swiss Tropical and Public Health Institute

This code is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.

This script is used to make a CSV file of intervention values

Author: Aatreyee Mimi Das
"""

#%%
# ----------- Import modules needed -------------
from itertools import product
import numpy as np
import pandas as pd
import os

#%%
# ------------- Set arrays of values for the parameters -------------------
follow_up = [0.0, 0.35, 1.0] # proportion of malaria cases followed up at the household level by DMSOs
neighbours = [0, 20, 100] # When following up at the household level, how many neighbours do we follow up? (0=None, household only)
treatment_seeking = [1, 2, 3] # Factor by which to multiply the current treatment seeking rate (1 = keep as is, no intervention)
rfmda = [0, 1] # Do we switch from RCD to rfMDA (change the test sensitivity to 1)? 0=No, 1=Yes
imp_treat = [0.0, 0.25, 0.5, 0.75, 0.90, 1.0] # proportion of imported cases that we treat early enough that they cannot lead to further cases
beta = [0.0, 0.5, 0.75, 0.9] # proportion by which the vectorial capacity, and thus beta, is reduced

tot_length = len(follow_up)*len(neighbours)*len(treatment_seeking)*len(rfmda)*len(imp_treat)*len(beta)
#tot_length = len(follow_up)*len(neighbours)*len(treatment_seeking)*len(rfmda)*len(imp_treat)
interventions_tb = np.zeros(shape=(tot_length,6)) #change to 6 when including beta

number = 0
for combination in product(follow_up, neighbours, treatment_seeking,\
                          rfmda, imp_treat, beta):

    interventions_tb[number,:] = combination
    number += 1

#%%
# ------------- Create full factorial table of combinations -------------------
# ------------------ and remove ones that are redundant -----------------------  
    
int_full_fac_df = pd.DataFrame(interventions_tb)
int_full_fac_df.columns = ['followup','neighbours','treatmentseeking','rfmda','imptreat','betared']

int_df = int_full_fac_df.loc[(int_full_fac_df['followup'] != 0.0) | \
                              (int_full_fac_df['rfmda'] != 1)]
int_df = int_df.loc[(int_df['followup'] != 0.0) | \
                              (int_df['neighbours'] != 20)]
int_df = int_df.loc[(int_df['followup'] != 0.0) | \
                              (int_df['neighbours'] != 100)]
int_df = int_df.iloc[np.repeat(np.arange(len(int_df)), 2)] #uncomment when not doing parameter uncertainty runs


os.chdir('~/interventions/')
int_df.to_csv('interventions.csv', index=False)
