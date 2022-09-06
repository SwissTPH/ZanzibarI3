#u -*- coding: utf-8 -*-
"""
Copyright (C) 2022 Swiss Tropical and Public Health Institute

This malaria model is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.


Stocastic implementation of a compartmental model for 
imported, introduced and indigenous malaria cases in the presence of 
RCD and human mobility

This script is the base script that runs the majority of interventions.

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
imptreat_separate = True

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

if imptreat_separate:
    df_int = pd.read_csv('interventions_imptreat.csv')
else:
    df_int = pd.read_csv('interventions.csv')

#%%
# ---------- Set up the simulations and the fixed parameter values ------------

n = 3 # Pemba, Unguja, Mainland Tanzania (this is the order used in later arrays and matrices)
nDays = 1  # frequency of RCD implementation

if server:
    tot_sim = 50. # total number of years over which to simulate
    init_sim = 10. # initial number of years over which the model is run to allow the number of indigenous, introduced, and imported cases to reach equilibrium
else:
    tot_sim = 10. # total number of years over which to simulate
    init_sim = 5. # initial number of years over which the model is run to allow the number of indigenous, introduced, and imported cases to reach equilibrium

if server:
    no_runs = 250 # to allow for running within 30 mins
else:
    no_runs = 1

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

pcr_pos_index_HH = np.asarray(df['pcr_pos_index_HH']) # prop. PCR-positive in the index HH, this does not include the index case
pcr_pos_neighbour = np.asarray(df['pcr_pos_neighbour']) # prop. PCR-positive in the neighbouring HHs
rdt_sensitivity = np.asarray(df['rdt_sensitivity']) # for RCD, this is 0.34, for rfMDA, this changes later on to 1
nu_HH = np.asarray(df['nu_HH']) # number of people investigated within index HH (including index case)
nu_n = np.asarray(df['nu_n']) # number of people investigated in neighbouring houses

pcr_pos_index_HH = (pcr_pos_index_HH*(nu_HH-1) + 1)/nu_HH # prevalence including the index case - to ensure that index cases are treated
tau_HH = pcr_pos_index_HH/I_eq  # targetting ratio within index households
tau_n = pcr_pos_neighbour/I_eq  # targetting ratio in neighbours who are also followed up

cases_per_day = np.asarray(df['cases_per_day']) # cases arriving at any health facility per day on each island
treatment_seeking = cases_per_day/(I_eq*pop) # rate at which an infected person seeks treatment per day
followup_days = 3 # the number of days in which the cases are followed up (determines the prop. of cases followed up), can be 3, 6, 15, or 21

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
# ------------ For an array job, take the run number for the bash -------------
# ------------ script and use this to determine the interventions -------------
if server:
    run_number = int(argv[1])-1
else:
    run_number = 1

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

beta = np.array(beta)[0]   # This is needed to get a 1D array

#%%
# ---------- This function is the bulk of the Gillespie algorithm -------------
# ------ It calculates the number of state transitions in each timestep -------
# ---------- based on a binomial tau-leap simulation algorithm ----------------

@njit
def stoc_eqs(INPUT, lop, cumsum, cumrec, beta, treatment_seeking, followup_prop, phi_consts,
             nu, theta_in, theta_out): 
    Rate = np.zeros((9*n))
    Change = np.zeros((6*n,4))
    Change = Change.astype(np.int64)
    S = np.copy(INPUT[0,:])
    P = np.copy(INPUT[1,:])
    T = np.copy(INPUT[2,:])
    D = np.copy(INPUT[3,:])
    cP = np.copy(cumsum[0:n])
    cT = np.copy(cumsum[n:2*n])
    cD = np.copy(cumsum[2*n:3*n])
    rP = np.copy(cumrec[0:n])
    rT = np.copy(cumrec[n:2*n])
    rD = np.copy(cumrec[2*n:3*n])
    rdts = np.copy(cumrec[3*n:4*n])
    acts = np.copy(cumrec[4*n:5*n])
    
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
        # Natural clearance
        Rate[3*n+k] = mu[k]*P[k]; Change[3*n+k,:] = [+1,-1,0,0]
        Rate[4*n+k] = mu[k]*T[k]; Change[4*n+k,:] = [+1,0,-1,0]
        Rate[5*n+k] = mu[k]*D[k]; Change[5*n+k,:] = [+1,0,0,-1]
        
        # RCD 
        Rate[6*n+k] = nDays*phi[k]*P[k]/N[k]
        Rate[7*n+k] = nDays*phi[k]*T[k]/N[k]
        Rate[8*n+k] = nDays*phi[k]*D[k]/N[k]        
        

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
        if P[j]<0 or T[j]<0 or D[j]<0:
            print('Error, below 0',P,T,D)
        # imported recovery events
        if P[j]==0:
            p1 = 0
            p2 = 0
        else:
            p1 = Rate[3*n+j]/P[j]
            p2 = Rate[6*n+j]/P[j]
        if p1>1:
            p1 = 1
        if p2>1:
            p2 = 1
        Draws3 = np.random.binomial(P[j],p1+p2,1)  # draw for all recovery events (natural clearance + RCD)
        if p1+p2==0 or Draws3[0]==0:
            Draws6 = np.random.binomial(1,0,1)
        else:
            Draws6 = np.random.binomial(Draws3[0],p2/(p1+p2),1)  # draw again for RCD recovery events
        INPUT[:,j] += Draws3*Change[3*n+j,:]
        rP[j] += Draws3[0]
        acts[j] += Draws6[0]
        
        # introduced recovery events
        if T[j]==0:
            p1 = 0
            p2 = 0
        else:
            p1 = Rate[4*n+j]/T[j]
            p2 = Rate[7*n+j]/T[j]
        if p1>1:
            p1=1
        if p2>1:
            p2=1
        Draws4 = np.random.binomial(T[j],p1+p2,1)  # draw for all recovery events (natural clearance + RCD)
        if p1+p2==0 or Draws4[0]==0:
            Draws7 = np.random.binomial(1,0,1)
        else:
            Draws7 = np.random.binomial(Draws4[0],p2/(p1+p2),1)  # draw again for RCD recovery events
        INPUT[:,j] += Draws4*Change[4*n+j,:]
        rT[j] += Draws4[0]
        acts[j] += Draws7[0]
        
        # indigenous recovery events
        if D[j]==0:
            p1 = 0
            p2 = 0
        else:
            p1 = Rate[5*n+j]/D[j]
            p2 = Rate[8*n+j]/D[j]
        if p1>1:
            p1=1
        if p2>1:
            p2=1
        Draws5 = np.random.binomial(D[j],p1+p2,1)
        if p1+p2==0 or Draws5[0]==0:
            Draws8 = np.random.binomial(1,0,1)
        else:
            Draws8 = np.random.binomial(Draws5[0],p2/(p1+p2),1)           
        INPUT[:,j] += Draws5*Change[5*n+j,:]
        
        rD[j] += Draws5[0]
        acts[j] += Draws8[0]
        
        rdts[j] += iota[j]*nu[j]  # number of people followed up x number of people tested per case followed up
            
    cumsum = np.concatenate((cP,cT,cD))
    cumrec = np.concatenate((rP,rT,rD,rdts,acts))
    return INPUT, cumsum, cumrec

#%%
# ------- This function runs the simulation repeatedly for the whole time -----
# ----------------- period, and updates the outputs arrays --------------------

# First we run for 10 years to allow the model to reach equilibrium, as 
# the initial conditions put all cases as indigenous cases.
# Then, after 10 years, we introduce interventions.
    
@njit
def Stoch_Iteration(INPUT, cumsum, cumrec):
    # State variables
    S = np.zeros(shape=(len(Time_year),n))
    P = np.zeros(shape=(len(Time_year),n))
    T = np.zeros(shape=(len(Time_year),n))
    D = np.zeros(shape=(len(Time_year),n))
    # Cumulative incidence
    cP = np.zeros(shape=(len(Time_year),n))
    cT = np.zeros(shape=(len(Time_year),n))
    cD = np.zeros(shape=(len(Time_year),n))
    # Cumulative recovery
    rP = np.zeros(shape=(len(Time_year),n))
    rT = np.zeros(shape=(len(Time_year),n))
    rD = np.zeros(shape=(len(Time_year),n))
    # Tracking RDT and ACT use
    rdts = np.zeros(shape=(len(Time_year),n))
    acts = np.zeros(shape=(len(Time_year),n))
    
    i = 0
    for lop in Time:
        if lop < 365*init_sim: # if time is within the initial time needed to reach equilibrium, run without interventions
            res,cs,cr = stoc_eqs(INPUT, lop, cumsum, cumrec, beta, treatment_seeking, followup_prop, phi_consts,
             nu=nu_HH, theta_in=theta, theta_out=theta)
            INPUT = res
            cumsum = cs
            cumrec = cr
            
            for ss in Time_year:  # these two lines are included because 'in' for int in an array doesn't work with numba
                if lop==ss:
                    S[i,:] = INPUT[0,:]
                    P[i,:] = INPUT[1,:]
                    T[i,:] = INPUT[2,:]
                    D[i,:] = INPUT[3,:]
                    cP[i,:] = cumsum[0:n]
                    cT[i,:] = cumsum[n:2*n]
                    cD[i,:] = cumsum[2*n:3*n]
                    rP[i,:] = cumrec[0:n]
                    rT[i,:] = cumrec[n:2*n]
                    rD[i,:] = cumrec[2*n:3*n]
                    rdts[i,:] = cumrec[3*n:4*n] 
                    acts[i,:] = cumrec[4*n:5*n]
                    i += 1
            
        else: # else, run with interventions
            res,cs,cr = stoc_eqs(INPUT, lop, cumsum, cumrec, beta=beta_new, treatment_seeking=treatment_seeking_new, 
                              followup_prop=followup_new, phi_consts=phi_consts_new, nu=nu_new, 
                                 theta_in=theta_in, theta_out=theta_out)
            INPUT = res
            cumsum = cs
            cumrec = cr
            
            for ss in Time_year:
                if lop==ss:
                    S[i,:] = INPUT[0,:]
                    P[i,:] = INPUT[1,:]
                    T[i,:] = INPUT[2,:]
                    D[i,:] = INPUT[3,:]
                    cP[i,:] = cumsum[0:n]
                    cT[i,:] = cumsum[n:2*n]
                    cD[i,:] = cumsum[2*n:3*n]
                    rP[i,:] = cumrec[0:n]
                    rT[i,:] = cumrec[n:2*n]
                    rD[i,:] = cumrec[2*n:3*n]
                    rdts[i,:] = cumrec[3*n:4*n] 
                    acts[i,:] = cumrec[4*n:5*n]
                    i += 1
            
    return S,P,T,D,cP,cT,cD,rP,rT,rD,rdts,acts


#%%    
# ----------- Interventions ----------------

followup_new = df_int['followup'][run_number]
neighbours_new = df_int['neighbours'][run_number]
treatmentseeking_factor = df_int['treatmentseeking'][run_number]
rfmda_new = df_int['rfmda'][run_number]
beta_red = df_int['betared'][run_number]

if imptreat_separate:
    imptreat_cp_new = df_int['imptreat_cp'][run_number]
    imptreat_toa_new = df_int['imptreat_toa'][run_number]
else:
    imptreat_new = df_int['imptreat'][run_number]


tau_times_nu = tau_HH*nu_HH+tau_n*neighbours_new # no. of neighbours range from 0 to 100
nu_new = nu_HH+neighbours_new  # measuring this to keep track of how many RDTs are used in RCD

# treatment seeking
treatment_seeking_new = treatment_seeking*treatmentseeking_factor 

# rfMDA
if rfmda_new == 1:
    rdt_sensitivity_new = 1 # if TRUE, then no need to test before treating, so no cases missed due to low RDT sensitivity
else:
    rdt_sensitivity_new = rdt_sensitivity

phi_consts_new = tau_times_nu*rdt_sensitivity_new

# importation treatment
if imptreat_separate:
    A = imptreat_toa_new
    B = imptreat_cp_new
else:
    A = imptreat_new
    B = imptreat_new

theta_out = theta.copy() #.copy()used because otherwise a python assigns a pointer, which changes theta as well
theta_out[0:2,2] = (1-A)*theta[0:2,2]
theta_in = theta.copy()
theta_in[2,0:2] = (1-B)*theta[2,0:2]

# reduction in vectorial capacity
beta_new = beta*np.array([1-beta_red, 1-beta_red, 1])

#%%
# ------------------ Set up initial values -----------------------------------

P0 = np.zeros((n)) # initial imported cases (assuming none)
T0 = np.zeros((n)) # initial introduced cases (assuming none)
D0 = np.round(pop*I_eq) # initial indigenous cases (assuming all current cases are indigenous cases)
S0 = pop-P0-T0-D0

# variables for cumulative cases and recoveries
cumsum = np.zeros((3*n))
cumrec = np.zeros((5*n))

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
rP_out = np.zeros(shape = (len(Time_year)*no_runs,n))
rT_out = np.zeros(shape = (len(Time_year)*no_runs,n))
rD_out = np.zeros(shape = (len(Time_year)*no_runs,n))
rdts_out = np.zeros(shape = (len(Time_year)*no_runs,n))
acts_out = np.zeros(shape = (len(Time_year)*no_runs,n))

#%%
# ------------- Run simulations ------------------------

start = perf_counter()

seed_start = 0
if run_number%2==1: seed_start=no_runs # to have seeds running from 1 to 500, rather than 1 to 250 twice
np.random.seed(seed_start)

for k in range(no_runs):
    INPUT = np.array((S0,P0,T0,D0))
    S,P,T,D,cP,cT,cD,rP,rT,rD,rdts,acts = Stoch_Iteration(INPUT, cumsum, cumrec)
    P_out[k*len(Time_year):(k+1)*len(Time_year),:] = P
    T_out[k*len(Time_year):(k+1)*len(Time_year),:] = T
    D_out[k*len(Time_year):(k+1)*len(Time_year),:] = D
    I_out[k*len(Time_year):(k+1)*len(Time_year),:] = P+T+D
    cP_out[k*len(Time_year):(k+1)*len(Time_year),:] = cP
    cT_out[k*len(Time_year):(k+1)*len(Time_year),:] = cT
    cD_out[k*len(Time_year):(k+1)*len(Time_year),:] = cD
    rP_out[k*len(Time_year):(k+1)*len(Time_year),:] = rP
    rT_out[k*len(Time_year):(k+1)*len(Time_year),:] = rT
    rD_out[k*len(Time_year):(k+1)*len(Time_year),:] = rD
    rdts_out[k*len(Time_year):(k+1)*len(Time_year),:] = rdts
    acts_out[k*len(Time_year):(k+1)*len(Time_year),:] = acts
    
files_list = ['P_out_', 'T_out_', 'D_out_', 'I_out_', 'cP_out_', 'cT_out_', 'cD_out_',
              'rP_out_', 'rT_out_', 'rD_out_', 'rdts_out_', 'acts_out_']
outputs_list = [P_out, T_out, D_out, I_out, cP_out, cT_out, cD_out, rP_out, rT_out, 
                rD_out, rdts_out, acts_out]

#%%
# -------------- Save outputs -------------------
if imptreat_separate:
    chdir('../outputs/imptreat')
    for j in range(len(files_list)):
        filename = files_list[j] + str(run_number) + '.csv'
        with open(filename, 'w') as f:
            np.savetxt(f, outputs_list[j], fmt='%5s', delimiter=',')
        #f.close()
else:
    chdir('../outputs/i3_model')
    for j in range(len(files_list)):
        filename = files_list[j] + str(run_number) + '.csv'
        with open(filename, 'w') as f:
            np.savetxt(f, outputs_list[j], fmt='%5s', delimiter=',')
        #f.close()


end = perf_counter()
print("Time taken (seconds):", end-start)

