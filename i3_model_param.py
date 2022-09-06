# -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Swiss Tropical and Public Health Institute

This malaria model is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.


Stocastic implementation of a compartmental model for 
imported, introduced and indigenous malaria cases in the presence of 
RCD and human mobility

This script uses parameters draw from distributions to explore the parameter 
uncertainty in the outputs

Author: Aatreyee Mimi Das
"""
#%%
# ----------- Import modules needed -------------
import numpy as np
from numba import njit
import pandas as pd
from os import chdir
from time import perf_counter
from sys import argv

#%%
# ------------- Load data and intervention parameters -------------------
server = True 
# False option is for quickly testing that the model is running on your local machine

# Directories need to be adapted to where files are stored
if server:
    chdir('~/data/')
else:    
    chdir('~/data/')

df = pd.read_csv('RADZEC_data.csv')


if server:
    chdir('~/interventions/')
else:    
    chdir('~/interventions/')

df_int = pd.read_csv('interventions_param.csv')
    
#%% 
# ------------ For an array job, take the run number for the bash -------------
# ------------ script and use this to determine the interventions -------------
if server:
    run_number = int(argv[1])-1
else:
    run_number = 15731

#%%
# ------------- Load parameter draws from distributions  -------------------
df_params = pd.read_csv('param_100.csv')
df_params = df_params.drop(['Unnamed: 0'], axis=1)
parcols = df_params.columns.values

par = df_params.loc[run_number % 100]
par = par.values    
print(par)

#%%
# ---------- Set up the simulations and the fixed parameter values ------------

n = 3 # Pemba, Unguja, Mainland Tanzania (this is the order used in later arrays and matrices)
nDays = 1  # frequency of RCD implementation
tot_sim = 50. # total number of years over which to simulate
init_sim = 10. # initial number of years over which the model is run to allow the number of indigenous, introduced, and imported cases to reach equilibrium
no_runs = 5 # no. runs for each parameter set

mu = np.asarray(df['mu'])  # currently assuming 200 days for infection clearance (no treatment)
pop = np.asarray(df['pop'])  # Total population in each patch

# Decide whether movement from the mainland will be included or not - for all figures presented, this was set to 'True'
movement_mainland = True
if movement_mainland:
    theta = np.mat('0.990698 0.003869 0.000057; 0.003182 0.970198 0.000533; 0.00612 0.025933 0.99941')
else:
    theta = np.mat('0.990698 0.003869 0; 0.003182 0.970198 0; 0.00612 0.025933 1')

#%%
# ---- Decide whether index cases or index households or neither -------
# ---------- are to be included in the baseline prevalence -------------

index_cases = False  # everyone included including index cases and households
index_house = False  # index cases not included but index households are included
neither = True  # neither index cases nor households are included

if index_cases:
    I_eq = df['I_eq_index']
elif index_house:
    I_eq = df['I_eq_index_HH']
elif neither:
    I_eq = df['I_eq_neither']
else:
    print('Error: one of the three options needs to be selected.')

I_eq = np.asarray(I_eq)

#%%
# ---------- Set the RCD parameters -----------------

rdt_sensitivity = np.asarray(df['rdt_sensitivity']) # For RCD, this is 0.34, for rfMDA, this changes later on to 1

# The following parameters have been drawn from distributions
I_eq[0] = par[0]
I_eq[1] = par[1]
nu_HH = np.hstack((par[6:8],[1]))
nu_n = np.hstack((par[8:10],[1]))
pcr_pos_index_HH = np.hstack((par[2:4],[0])) 
pcr_pos_neighbour = np.hstack((par[4:6],[0]))

pcr_pos_index_HH = (pcr_pos_index_HH*(nu_HH-1) + 1)/nu_HH # prevalence including the index case - to ensure that index cases are treated
tau_HH = pcr_pos_index_HH/I_eq  # targetting ratio within index households
tau_n = pcr_pos_neighbour/I_eq  # targetting ratio in neighbours who are also followed up

cases_per_day = np.asarray(df['cases_per_day']) # Cases arriving at any health facility per day on each island
treatment_seeking = cases_per_day/(I_eq*pop) # rate at which an infected person seeks treatment per day
followup_days = 3 # number of days of follow up allowed can be 3, 6, 15, or 21

if followup_days==3:  # change followup_prop to proportion
    followup_prop = np.asarray(df['followup_3_days'])
elif followup_days==6:
    followup_prop = np.asarray(df['followup_6_days'])
elif followup_days==15:
    followup_prop = np.asarray(df['followup_15_days'])
elif followup_days==21:
    followup_prop = np.asarray(df['followup_21_days'])
else:
    print('Error: number of days of follow up must be 3, 6, 15 or 21.')

phi_consts = rdt_sensitivity*(tau_HH*nu_HH)

#%%   
# -------- Back calculation of the transmission parameter beta ----------------
    
I_eff = np.zeros((n))
for i in range(n):
    I_eff[i] = np.dot(pop*I_eq, np.transpose(theta[i,:]))/np.dot(pop,np.transpose(theta[i,:]))
A = np.diag(I_eff)
M = np.transpose(np.dot(A,theta))

iota = I_eq*pop*treatment_seeking*followup_prop
phi = phi_consts*iota*I_eq/pop
phi[2] = 0  # RCD not present on the mainland

p = (mu*I_eq + phi)/(1-I_eq)

# Calculate beta using the inverse of M and p

M_inv = np.linalg.inv(M)
beta = np.dot(M_inv,np.transpose(p))
print('Beta for Pemba, Unguja and Mainland Tanzania\n',beta)
beta = np.array(beta)[0]   

#%%
# ---------- This function is the bulk of the Gillespie algorithm -------------
# ------ It calculates the number of state transitions in each timestep -------
# ---------- based on a binomial tau-leap simulation algorithm ----------------

@njit
def stoc_eqs(INPUT, lop, cumsum, beta, treatment_seeking, followup_prop, phi_consts,
             theta_in, theta_out): 
    Rate = np.zeros((6*n))
    Change = np.zeros((6*n,4))
    Change = Change.astype(np.int64)
    S = np.copy(INPUT[0,:])
    P = np.copy(INPUT[1,:])
    T = np.copy(INPUT[2,:])
    D = np.copy(INPUT[3,:])
    cP = np.copy(cumsum[0:n])
    cT = np.copy(cumsum[n:2*n])
    cD = np.copy(cumsum[2*n:3*n])
    
    I_ef = np.zeros((n)) # could not take I_eff from before because it was a global variable and numba couldn't compile it, so new variable made
    I_ef2 = np.zeros((n))

    I = P+T+D # total number infected in absolute terms

    N = S+I


        
    if lop%nDays==0:
        iota = I*treatment_seeking*followup_prop
        phi = phi_consts*iota
    else:
        phi = np.zeros((n))
    phi[2] = 0 # ensure no RCD on mainland
    

    for k in range(n):
        # TRANSMISSION EVENTS
        # importation transmission events
        
        for i in range(n):
            if i==k:
                I_ef[i] = 0
            else:
                I_ef[i]=np.dot(I, np.transpose(theta_out[i,:]))/np.dot(N,np.transpose(theta[i,:]))
        Rate[k] = np.multiply(S[k], np.dot(beta*I_ef, theta_in[:,k])); Change[k,:] = [-1,+1,0,0]
        
        # introduction transmission events
        I_ef2[k]=(np.dot(I, np.transpose(theta_out[k,:]))-I[k]*theta_out[k,k])/np.dot(N,np.transpose(theta[k,:]))
        Rate[n+k] = beta[k]*(P[k]*theta_out[k,k]/np.dot(N,np.transpose(theta[k,:]))+I_ef2[k])*theta_in[k,k]*S[k]
        Change[n+k,:] = [-1,0,+1,0]
        
        # indigenous transmission events
        Rate[2*n+k] = beta[k]*(theta_out[k,k]/np.dot(N,np.transpose(theta[k,:])))*(T[k]+D[k])*theta_in[k,k]*S[k]
        Change[2*n+k,:] = [-1,0,0,+1]
        
        # RECOVERY EVENTS
        # Includes natural clearance and RCD
        Rate[3*n+k] = (mu[k]*N[k]+nDays*phi[k])*P[k]/N[k]; Change[3*n+k,:] = [+1,-1,0,0]
        Rate[4*n+k] = (mu[k]*N[k]+nDays*phi[k])*T[k]/N[k]; Change[4*n+k,:] = [+1,0,-1,0]
        Rate[5*n+k] = (mu[k]*N[k]+nDays*phi[k])*D[k]/N[k]; Change[5*n+k,:] = [+1,0,0,-1]

    for j in range(n):
        # importation transmission events
        if S[j]==0:
            p = 0
        else:
            p = Rate[j]*TS/S[j]
        if p>1:
            p=1
        Draws0 = np.random.binomial(S[j],p,1)
        INPUT[:,j] += Draws0*Change[j,:]
        cP[j] += Draws0[0]
        # introduction transmission events
        if S[j]==0:
            p = 0
        else:
            p = Rate[n+j]*TS/S[j]
        if p>1:
            p=1
        Draws1 = np.random.binomial(S[j],p,1)
        INPUT[:,j] += Draws1*Change[n+j,:]
        cT[j] += Draws1[0]
        # indigenous transmission events
        if S[j]==0:
            p = 0
        else:
            p = Rate[2*n+j]*TS/S[j]
        if p>1:
            p=1
        Draws2 = np.random.binomial(S[j],p,1)
        INPUT[:,j] += Draws2*Change[2*n+j,:]
        cD[j] += Draws2[0]
        
        
    for j in range(n):
        # imported recovery events
        if P[j]==0:
            p = 0
        else:
            p = Rate[3*n+j]/P[j]
        if p>1:
            p=1
        Draws3 = np.random.binomial(P[j],p,1)
        INPUT[:,j] += Draws3*Change[3*n+j,:]
        # introduced recovery events
        if T[j]==0:
            p = 0
        else:
            p = Rate[4*n+j]/T[j]
        if p>1:
            p=1
        Draws4 = np.random.binomial(T[j],p,1)
        INPUT[:,j] += Draws4*Change[4*n+j,:]
        # indigenous recovery events
        if D[j]==0:
            p = 0
        else:
            p = Rate[5*n+j]/D[j]
        if p>1:
            p=1
        Draws5 = np.random.binomial(D[j],p,1)
        INPUT[:,j] += Draws5*Change[5*n+j,:]
    cumsum = np.concatenate((cP,cT,cD))
    return INPUT, cumsum

#%%
# ------- This function runs the simulation repeatedly for the whole time -----
# ----------------- period, and updates the outputs arrays --------------------

# First we run for 10 years to allow the model to reach equilibrium, as 
# the initial conditions put all cases as indigenous cases.
# Then, after 10 years, we introduce interventions.    


@njit
def Stoch_Iteration(INPUT, cumsum):
    # State variables
    S = np.zeros(shape=(len(Time_year),n))
    P = np.zeros(shape=(len(Time_year),n))
    T = np.zeros(shape=(len(Time_year),n))
    D = np.zeros(shape=(len(Time_year),n))
    # Cumulative incidence
    cP = np.zeros(shape=(len(Time_year),n))
    cT = np.zeros(shape=(len(Time_year),n))
    cD = np.zeros(shape=(len(Time_year),n))
    i = 0
    for lop in Time:
        if lop < 365*init_sim: # if time is within the initial time needed to reach equilibrium, run without interventions
            res,cs = stoc_eqs(INPUT, lop, cumsum, beta, treatment_seeking, followup_prop, phi_consts,
             theta_in=theta, theta_out=theta)
            for ss in Time_year:  # these two lines are included because 'in' for int in an array doesn't work with numba
                if lop==ss:
                    S[i,:] = INPUT[0,:]
                    P[i,:] = INPUT[1,:]
                    T[i,:] = INPUT[2,:]
                    D[i,:] = INPUT[3,:]
                    cP[i,:] = cumsum[0:n]
                    cT[i,:] = cumsum[n:2*n]
                    cD[i,:] = cumsum[2*n:3*n]
                    i += 1
            INPUT = res
            cumsum = cs
        else: # else, run with interventions
            res,cs = stoc_eqs(INPUT, lop, cumsum, beta=beta, treatment_seeking=treatment_seeking_new, 
                              followup_prop=followup_new, phi_consts=phi_consts_new,
             theta_in=theta_in, theta_out=theta_out)
            for ss in Time_year:
                if lop==ss:
                    S[i,:] = INPUT[0,:]
                    P[i,:] = INPUT[1,:]
                    T[i,:] = INPUT[2,:]
                    D[i,:] = INPUT[3,:]
                    cP[i,:] = cumsum[0:n]
                    cT[i,:] = cumsum[n:2*n]
                    cD[i,:] = cumsum[2*n:3*n]
                    i += 1
            INPUT = res
            cumsum = cs
    return S,P,T,D,cP,cT,cD


#%%    
# ----------- Interventions ----------------

followup_new = df_int['followup'][run_number // 100]
neighbours_new = df_int['neighbours'][run_number // 100]
treatmentseeking_factor = df_int['treatmentseeking'][run_number // 100]
rfmda_new = df_int['rfmda'][run_number // 100]
imptreat_new = df_int['imptreat'][run_number // 100]

tau_times_nu = tau_HH*nu_HH+tau_n*neighbours_new # no. of neighbours range from 0 to 100

treatment_seeking_new = treatment_seeking*treatmentseeking_factor 
    
if rfmda_new == 1:
    rdt_sensitivity_new = 1 # if TRUE, then no need to test before treating, so no cases missed due to low RDT sensitivity
else:
    rdt_sensitivity_new = rdt_sensitivity
    
phi_consts_new = tau_times_nu*rdt_sensitivity_new

# importation treatment
A = imptreat_new
B = imptreat_new

theta_out = theta.copy() #.copy()used because otherwise a python assigns a pointer, which changes theta as well
theta_out[0:2,2] = (1-A)*theta[0:2,2]
theta_in = theta.copy()
theta_in[2,0:2] = (1-B)*theta[2,0:2]

#%%
# ------------------ Set up initial values -----------------------------------

P0 = np.zeros((n)) # initial imported cases (assuming none)
T0 = np.zeros((n)) # initial introduced cases (assuming none)
D0 = np.round(pop*I_eq) # initial indigenous cases (assuming all current cases are indigenous cases)
S0 = pop-P0-T0-D0

# variables for cumulative cases
cumsum = np.zeros((3*n))

# Run simulations
MaxTime = 365*tot_sim
TS = 1 # timestep tau for tau-leap
Time=np.arange(0, MaxTime, TS)
Time_year=np.arange(365, MaxTime+365, 365/TS)-1

# Set up arrays for recording outputs
P_out = np.zeros(shape = (len(Time_year)*no_runs,n))
T_out = np.zeros(shape = (len(Time_year)*no_runs,n))
D_out = np.zeros(shape = (len(Time_year)*no_runs,n)) 
I_out = np.zeros(shape = (len(Time_year)*no_runs,n))
cP_out = np.zeros(shape = (len(Time_year)*no_runs,n))
cT_out = np.zeros(shape = (len(Time_year)*no_runs,n))
cD_out = np.zeros(shape = (len(Time_year)*no_runs,n))

#%%
# ------------- Run simulations ------------------------

start = perf_counter()

seed_start = 0
np.random.seed(seed_start)


for k in range(no_runs):
    INPUT = np.array((S0,P0,T0,D0))
    S,P,T,D,cP,cT,cD = Stoch_Iteration(INPUT, cumsum)
    P_out[k*len(Time_year):(k+1)*len(Time_year),:] = P
    T_out[k*len(Time_year):(k+1)*len(Time_year),:] = T
    D_out[k*len(Time_year):(k+1)*len(Time_year),:] = D
    I_out[k*len(Time_year):(k+1)*len(Time_year),:] = P+T+D
    cP_out[k*len(Time_year):(k+1)*len(Time_year),:] = cP
    cT_out[k*len(Time_year):(k+1)*len(Time_year),:] = cT
    cD_out[k*len(Time_year):(k+1)*len(Time_year),:] = cD
    
files_list = ['P_out_', 'T_out_', 'D_out_', 'I_out_', 'cP_out_', 'cT_out_', 'cD_out_']
outputs_list = [P_out, T_out, D_out, I_out, cP_out, cT_out, cD_out]

#%%
# -------------- Save outputs -------------------
if server:
    chdir('../outputs/params')
    for j in range(len(files_list)):
        filename = files_list[j] + str(run_number) + '.csv'
        print(filename)
        with open(filename, 'w') as f:
            np.savetxt(f, outputs_list[j], fmt='%5s', delimiter=',')
        f.close()


end = perf_counter()
print("Time taken (seconds):", end-start)
