# u -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Swiss Tropical and Public Health Institute

This malaria model is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.


Stocastic implementation of a compartmental model for 
imported, introduced and indigenous malaria cases in the presence of 
RCD and human mobility

This script runs the model with parameter values drawn from a Latin Hypercube in order to calculate Sobol indices

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
    chdir('~/sensitivity_analysis/')
else:    
    chdir('~/sensitivity_analysis/')

# read in sensitivity analysis parameter values
n_samples = 32768
sample_type = 'saltelli'  # can be 'LHS' or 'saltelli'. Note, only 'saltelli' works with SALib for Sobol indices

filename = 'param_values_'+ sample_type +'_' + str(n_samples) + '.txt'

df_sa_full = pd.DataFrame(np.loadtxt(filename))
df_sa_full.columns = ['betaP', 
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
              'thetaUM']

#%%
# ------------ For an array job, take the run number for the bash -------------
# ------------ script and use this to determine the interventions -------------
if server:
    run_number = int(argv[1])-1
else:
    run_number = 1

if server:
    no_runs = 128 # 128 runs per job so we can run in the 30 min queue
else:
    no_runs = 10 #smaller for testing

sample_no = run_number*no_runs

     
#%%
# ---------- Set up the simulations and the fixed parameter values ------------

n = 3 # Pemba, Unguja, Mainland Tanzania (this is the order used in later arrays and matrices)
nDays = 1  # frequency of RCD implementation
tot_sim = 50. # total number of years over which to simulate


pop = np.asarray(df['pop'])  # Total population in each patch



#%%
# ---------- Set the RCD parameters -----------------

nu_HH = np.asarray(df['nu_HH']) # number of people investigated within index HH (including index case)
nu_n = np.asarray(df['nu_n']) # number of people investigated in neighbouring houses


#%%
# ------ This function is for the sensitivity anaysis parameters---------------

def sa_params(df_sa_full, k):
    df_sa = df_sa_full.iloc[[k]]
    
    beta = np.empty(n)
    inflength = np.empty(n)
    tau_HH = np.empty(n)
    treatment_seeking = np.empty(n)
    rdt_sensitivity = np.empty(n)
    followup_prop = np.empty(n)
    theta = np.zeros((n,n))

    # Transmission and recovery rates
    beta[0] = df_sa['betaP']
    beta[1] = df_sa['betaU']
    beta[2] = df_sa['betaM'] 

    inflength.fill(df_sa['inflength'].values[0]) # the [0] is to select the only element so it can be used with .fill()
    mu  = 1/inflength
    
    # RCD parameters
    tau_HH = np.array([df_sa['tauHHP'].values[0], df_sa['tauHHU'].values[0], 0])
    tau_HH_adjusted = (tau_HH*(nu_HH-1)+1)/nu_HH
    treatment_seeking = np.array([df_sa['treatseekP'].values[0], df_sa['treatseekU'].values[0], 0])
    rdt_sensitivity.fill(df_sa['rdtsens'].values[0])
    followup_prop.fill(df_sa['followup'].values[0])
    phi_consts = rdt_sensitivity*(tau_HH_adjusted*nu_HH)
    
    # The movement matrix is an n x n matrix where all of the columns sum to 1
    theta[0,1] = df_sa['thetaPU']
    theta[0,2] = df_sa['thetaPM']
    theta[1,0] = df_sa['thetaUP']
    theta[1,2] = df_sa['thetaUM']
    theta[2,0] = df_sa['thetaMP']
    theta[2,1] = df_sa['thetaMU']
    
    # Time spent at home is determined by time spent away
    theta[0,0] = 1 - theta[1,0] - theta[2,0]
    theta[1,1] = 1 - theta[0,1] - theta[2,1]
    theta[2,2] = 1 - theta[0,2] - theta[1,2]
    
    
    
    return beta, mu, treatment_seeking, rdt_sensitivity, followup_prop, phi_consts, theta

#%%
# ---------- This function is the bulk of the Gillespie algorithm -------------
# ------ It calculates the number of state transitions in each timestep -------

@njit
def stoc_eqs(INPUT, lop, cumsum, beta, mu, treatment_seeking, followup_prop, phi_consts,
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


        
    if lop%nDays==0: # on every nDays, RCD is implemented
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

@njit
def Stoch_Iteration(INPUT, cumsum, beta, mu, treatment_seeking, rdt_sensitivity, 
                    followup_prop, phi_consts, theta):
    S = np.zeros((n))
    P = np.zeros((n))
    T = np.zeros((n))
    D = np.zeros((n))
    cP =np.zeros((n))
    cT = np.zeros((n))
    cDF = np.zeros((n))
    cDP = np.zeros((n))
    Di = np.zeros((n))

    for lop in Time:
        res,cs = stoc_eqs(INPUT, lop, cumsum, beta, mu, treatment_seeking, followup_prop, phi_consts,
         theta_in=theta, theta_out=theta)
        if lop==(Time_final-365):
            cDP = cumsum[2*n:3*n]
        if lop==Time_final:
            S = INPUT[0,:]
            P = INPUT[1,:]
            T = INPUT[2,:]
            D = INPUT[3,:]
            cP = cumsum[0:n]
            cT = cumsum[n:2*n]
            cDF = cumsum[2*n:3*n]
        INPUT = res
        cumsum = cs
    # incidence of indigenous infections in final year 
    # (cumulative final - cumulatve penultimate)
    Di = cDF - cDP 
    return S,P,T,D,cP,cT,cDF, Di


#%%
# ------------------ Set up initial values -----------------------------------

P0 = np.zeros((n)) # initial imported cases (assuming none)
T0 = np.zeros((n)) # initial introduced cases (assuming none)
D0 = np.round(pop*0.01) # initial indigenous cases (assuming 1% prevalence on Pemba and Unguja, 8% on mainland)
D0[2] = np.round(pop[2]*0.08)
S0 = pop-P0-T0-D0

# variables for cumulative cases
cumsum = np.zeros((3*n))

# Run simulations
MaxTime = 365*tot_sim
TS = 1 # timestep tau for tau-leap
Time=np.arange(0, MaxTime, TS)
Time_final=Time[-1]

# Set up arrays for recording outputs
P_out = np.zeros(shape = (no_runs,n))
T_out = np.zeros(shape = (no_runs,n))
D_out = np.zeros(shape = (no_runs,n)) 
I_out = np.zeros(shape = (no_runs,n))
cP_out = np.zeros(shape = (no_runs,n))
cT_out = np.zeros(shape = (no_runs,n))
cD_out = np.zeros(shape = (no_runs,n))
Di_out = np.zeros(shape = (no_runs,n))

#%%
# ------------- Run simulations ------------------------

start = perf_counter()

np.random.seed(sample_no)
for j in range(no_runs):
    k = sample_no + j
    print(k)
    beta, mu, treatment_seeking, rdt_sensitivity, followup_prop, \
    phi_consts, theta = sa_params(df_sa_full, k)

    INPUT = np.array((S0,P0,T0,D0))
    S,P,T,D,cP,cT,cD, Di = Stoch_Iteration(INPUT, cumsum, beta, mu, treatment_seeking, 
                                       rdt_sensitivity, followup_prop, phi_consts, theta)
    I = P+T+D
    P_out[j,:] = P
    T_out[j,:] = T
    D_out[j,:] = D
    I_out[j,:] = I/pop
    cP_out[j,:] = cP
    cT_out[j,:] = cT
    cD_out[j,:] = cD
    Di_out[j,:] = Di
    
    
files_list = ['P_out_', 'T_out_', 'D_out_', 'I_out_', 'cP_out_', 'cT_out_', \
              'cD_out_', 'Di_out_']
outputs_list = [P_out, T_out, D_out, I_out, cP_out, cT_out, cD_out, Di_out]


#%%
# -------------- Save outputs -------------------
if server:
    chdir('~/outputs/SA/')
    for j in range(len(files_list)):
        filename = files_list[j] + sample_type + '_' + str(n_samples) +\
        '_run_no_' + str(run_number) + '.txt'
        with open(filename, 'w') as f:
            np.savetxt(f, outputs_list[j])
    #f.close()

end = perf_counter()
print("Time taken (seconds) for ", sample_type, ":", end-start)

