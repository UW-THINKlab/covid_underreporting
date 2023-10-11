#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 21:27:20 2022

@author: gracejia
"""

import sys
sys.path.append('/Users/gracejia/Documents/A-UW//covid19 NSF Project/TLI controller')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import InterventionSIR as sir
import parameters as params
import optimize_interventions as oi

np.random.seed(8544)

#%%
t_sim_max = 360
name = "maintain-suppress-time"
correct_rate = 0.75   # 0.05, 0.1, 0.25,0.5,0.75, 1
case = 3   # no underreporting -1
            # fixed rate        3

outpath = '/Users/gracejia/Documents/A-UW//covid19 NSF Project/TLI controller/out/'
taus = params.taus_figure_interventions
covid_sir = sir.InterventionSIR(
    b_func = sir.Intervention(),
    R0 = params.R0_default,
    gamma = params.gamma_default,
    inits = params.inits_default)
covid_sir.reset()
covid_sir.outpath = outpath

null_time, null_result = covid_sir.integrate_null(t_sim_max)
null_states = np.hstack((null_time.reshape((10000,1)), null_result))
pd.DataFrame(null_states, columns = ["time","S","I","R"]).to_csv(outpath + "null_states.csv")

## set intervention strategy and scenario parameters
covid_sir.case = case
covid_sir.b_func.strategy = name
covid_sir.b_func.correct_rate = correct_rate

## iterate over intervention durations
taus = [14,28,56]
initsOfCovid_sir = covid_sir.inits
 
for i_tau, tau in enumerate(taus):
    covid_sir.inits = initsOfCovid_sir
    covid_sir.case = case
    covid_sir.b_func.tau = tau
    S_i_expected = 0
    print("optimizing strategy for {} "
          "with tau = {}".format(name, tau))
    if name == "maintain-suppress-time":
        S_i_expected, f = oi.calc_Sf_opt(
            covid_sir.R0,
            covid_sir.gamma * tau)
        I_i_expected = covid_sir.I_of_S(S_i_expected)
        covid_sir.b_func.S_i_expected = S_i_expected
        covid_sir.b_func.I_i_expected = I_i_expected
        covid_sir.b_func.S_i_observed = covid_sir.inits[0]
        covid_sir.b_func.I_i_observed = covid_sir.inits[1]
        covid_sir.b_func.S_i_corrected = covid_sir.inits[0]
        covid_sir.b_func.I_i_corrected = covid_sir.inits[1]                
        covid_sir.b_func.f = f

    t_i_opt = covid_sir.t_of_S(S_i_expected)[0]
    covid_sir.b_func.t_i = t_i_opt
    
    covid_sir.reset()
    
    if name == "maintain-suppress-time":
        covid_sir.state = np.append(initsOfCovid_sir, initsOfCovid_sir)
        
    covid_sir.integrate(t_sim_max)
    times = np.array(covid_sir.time_ts)
    b_of_t =  np.array([
    covid_sir.b_func(time,
                     covid_sir.R0 * covid_sir.gamma,
                     covid_sir.gamma,
                     S,
                     I)
    for time, S, I in zip(
            times,
            covid_sir.state_ts[:, 0],
            covid_sir.state_ts[:, 1])]).astype("float")

    file_name = 'b(t) for case {} tau = {} cc = {}'.format(covid_sir.case, tau, correct_rate)
    df = pd.DataFrame({"time" : times, "b(t)" : b_of_t})
    df.to_csv(outpath + file_name,index = False)
 
        
#%% 
# read b_t data
file_name_baseline = ['b(t) for case -1 tau = 14 cc = 1',
             'b(t) for case -1 tau = 28 cc = 1',
             'b(t) for case -1 tau = 56 cc = 1',
             ]
file_name_nc_1 = ['b(t) for case 3 tau = 14 cc = 1',
             'b(t) for case 3 tau = 28 cc = 1',
             'b(t) for case 3 tau = 56 cc = 1',
             ]
file_name_c_5 = ['b(t) for case 3 tau = 14 cc = 0.05',
             'b(t) for case 3 tau = 28 cc = 0.05',
             'b(t) for case 3 tau = 56 cc = 0.05',
             ]
file_name_c_10 = ['b(t) for case 3 tau = 14 cc = 0.1',
             'b(t) for case 3 tau = 28 cc = 0.1',
             'b(t) for case 3 tau = 56 cc = 0.1',
             ]

file_name_c_25 = ['b(t) for case 3 tau = 14 cc = 0.25',
             'b(t) for case 3 tau = 28 cc = 0.25',
             'b(t) for case 3 tau = 56 cc = 0.25',
             ]
file_name_c_50 = ['b(t) for case 3 tau = 14 cc = 0.5',
             'b(t) for case 3 tau = 28 cc = 0.5',
             'b(t) for case 3 tau = 56 cc = 0.5',
             ]

file_name_c_75 = ['b(t) for case 3 tau = 14 cc = 0.75',
             'b(t) for case 3 tau = 28 cc = 0.75',
             'b(t) for case 3 tau = 56 cc = 0.75',
             ]

b_t_baseline = []
for i in range(0,len(file_name_baseline)):
    df = pd.read_csv(outpath + file_name_baseline[i])
    col_name_1 = "time_{}".format(7*2**i)
    col_name_2 = "duration = {}".format(7*2**i)
    df.rename(columns = {"time":col_name_1,"b(t)":col_name_2})
    b_t_baseline.append(df)

b_t_nc_1 = []
for i in range(0,len(file_name_nc_1)):
    df = pd.read_csv(outpath + file_name_nc_1[i])
    b_t_nc_1.append(df)
    
b_t_c_5 = []
for i in range(0,len(file_name_c_5)):
    df = pd.read_csv(outpath + file_name_c_5[i])
    b_t_c_5.append(df)

b_t_c_10 = []
for i in range(0,len(file_name_c_10)):
    df = pd.read_csv(outpath + file_name_c_10[i])
    b_t_c_10.append(df)

b_t_c_25 = []
for i in range(0,len(file_name_c_25)):
    df = pd.read_csv(outpath + file_name_c_25[i])
    b_t_c_25.append(df)

b_t_c_50 = []
for i in range(0,len(file_name_c_50)):
    df = pd.read_csv(outpath + file_name_c_50[i])
    b_t_c_50.append(df)

b_t_c_75 = []
for i in range(0,len(file_name_c_75)):
    df = pd.read_csv(outpath + file_name_c_75[i])
    b_t_c_75.append(df)


#%%
# plot the all b(t) for duration 56 days
plt.figure(figsize = (10,7)) 
plt.plot(b_t_baseline[2]['time'],b_t_baseline[2]['b(t)'] * 3.92*0.05, label = "baseline: no underreport, no correction")  
plt.plot(b_t_nc_1[2]['time'],b_t_nc_1[2]['b(t)'] * 3.92*0.05, label = "NC: underreport, no correction")
# plt.plot(b_t_c_5[2]['time'],b_t_c_5[2]['b(t)'] * 3.92*0.05, label = "C: underreport, correction rate = 0.05")
plt.plot(b_t_c_10[2]['time'],b_t_c_10[2]['b(t)'] * 3.92*0.05, label = "C: underreport, correction rate = 0.10")
plt.plot(b_t_c_25[2]['time'],b_t_c_25[2]['b(t)'] * 3.92*0.05, label = "C: underreport, correction rate = 0.25")      
plt.plot(b_t_c_50[2]['time'],b_t_c_50[2]['b(t)'] * 3.92*0.05, label = "C: underreport, correction rate = 0.50")
plt.plot(b_t_c_75[2]['time'],b_t_c_75[2]['b(t)'] * 3.92*0.05, label = "C: underreport, correction rate = 0.75")

plt.xlabel("time") 
plt.ylabel("transmission rate")
plt.title("Effect of Underreporting and Correction on Transmission rate for policy duration 56 days")
plt.legend()
plt.show()

#%%
#read states data
state_name_baseline = ["states for case -1 tau 56 cc 1.csv"]
state_name_nc = ["states for case 3 tau 56 cc 1.csv"]
state_name_c_5 = ["states for case 3 tau 56 cc 0.05.csv"]
state_name_c_10 = ["states for case 3 tau 56 cc 0.1.csv"]
state_name_c_25 = ["states for case 3 tau 56 cc 0.25.csv"]
state_name_c_50 = ["states for case 3 tau 56 cc 0.5.csv"]
state_name_c_75 = ["states for case 3 tau 56 cc 0.75.csv"]

state_baseline = []
state_nc = []
state_c_5 = []
state_c_10 = []
state_c_25 = []
state_c_50 = []
state_c_75 = []

state_baseline.append(pd.read_csv(outpath + state_name_baseline[0]))
state_nc.append(pd.read_csv(outpath + state_name_nc[0]))
state_c_5.append(pd.read_csv(outpath + state_name_c_5[0]))
state_c_10.append(pd.read_csv(outpath + state_name_c_10[0]))
state_c_25.append(pd.read_csv(outpath + state_name_c_25[0]))
state_c_50.append(pd.read_csv(outpath + state_name_c_50[0]))
state_c_75.append(pd.read_csv(outpath + state_name_c_75[0]))

timeseries = np.linspace(0,360,10000)
plt.figure(figsize = (10,7))
plt.plot(timeseries,state_baseline[0]['I'], label = "baseline: no underreport, no correction")  
plt.plot(timeseries,state_nc[0]['I'], label = "NC: underreport, no correction")
# plt.plot(timeseries,state_c_5[0]['I'], label = "C: underreport, correction rate = 0.05")
plt.plot(timeseries,state_c_10[0]['I'], label = "C: underreport, correction rate = 0.10")
plt.plot(timeseries,state_c_25[0]['I'], label = "C: underreport, correction rate = 0.25")      
plt.plot(timeseries,state_c_50[0]['I'], label = "C: underreport, correction rate = 0.50")
plt.plot(timeseries,state_c_75[0]['I'], label = "C: underreport, correction rate = 0.75")
plt.xlabel("time") 
plt.ylabel("Infectious individuals")
plt.title("Effect of Underreporting and Correction on Infectious for policy duration 56 days")
plt.legend(loc = "best")
plt.show()
#%%
# only plot the policy duration
# x1 = b_t_baseline[2]['time'][1900:4200]
# y1 = b_t_baseline[2]['b(t)'][1900:4200]
# x2 = b_t_nc_1[2]['time'][1900:4200]
# y2 = b_t_nc_1[2]['b(t)'][1900:4200]
# x3 = b_t_c_75[2]['time'][1900:4200]
# y3 = b_t_c_75[2]['b(t)'][1900:4200]
# x4 = b_t_c_25[2]['time'][1900:4200]
# y4 = b_t_c_25[2]['b(t)'][1900:4200]
# plt.figure(figsize = (10,7)) 
# plt.plot(x1,y1* 3.92*0.05, label = "baseline")  
# plt.plot(x2,y2 * 3.92*0.05, label = "NC")
# plt.plot(x3,y3* 3.92*0.05, label = "C, rate = 0.75")
# plt.plot(x4,y4 * 3.92*0.05, label = "C, rate = 0.25")      
# plt.xlabel("time") 
# plt.ylabel("transmission rate")
# plt.title("Effect of Underreporting and Correction on Transmission rate for policy duration 56 days")
# plt.legend()
# plt.show()



#%%
# read states data
state_name_c_75 = ["states for case 3 tau 56 cc 0.75.csv"]
state_name_c_25 = ["states for case 3 tau 56 cc 0.25.csv"]
state_c_75 = []
state_c_75.append(pd.read_csv(outpath + state_name_c_75[0]))
state_c_25 = []
state_c_25.append(pd.read_csv(outpath + state_name_c_25[0]))
    
plt.figure(figsize = (10,7)) 
# plt.plot(b_t_baseline[2]['time'],state_c_75[0]["S_obs"], label = "75_S_obs")  
# plt.plot(b_t_baseline[2]['time'],state_c_75[0]["S"], label = "75_S")
plt.plot(b_t_baseline[2]['time'],state_c_25[0]["S_obs"], label = "25_S_obs")
plt.plot(b_t_baseline[2]['time'],state_c_25[0]["S"], label = "25_S")      
plt.xlabel("time") 
plt.ylabel("State (Susceptible)")
plt.title("Effect of Underreporting and Correction on Suspetible Compartment for policy duration 56 days")
plt.legend()
plt.show()

#%%
# plot I compartment data 
plt.figure(figsize = (10,7)) 
# plt.plot(b_t_baseline[2]['time'],state_c_75[0]["S_obs"], label = "75_S_obs")  
# plt.plot(b_t_baseline[2]['time'],state_c_75[0]["S"], label = "75_S")
plt.plot(b_t_baseline[2]['time'],state_c_25[0]["I_obs"], label = "25_I_obs")
plt.plot(b_t_baseline[2]['time'],state_c_25[0]["I"], label = "25_I")      
plt.xlabel("time") 
plt.ylabel("State (Infectious)")
plt.title("Effect of Underreporting and Correction on Suspetible Compartment for policy duration 56 days")
plt.legend()
plt.show()
#%%
# plot only the correct rate 15%
file_name_c_15 = ['b(t) for case 3 tau = 56 cc = 0.15']
b_c_15 = []
b_c_15.append(pd.read_csv(outpath + file_name_c_15[0]))
plt.figure(figsize = (10,7)) 
# plt.plot(b_t_baseline[2]['time'],state_c_75[0]["S_obs"], label = "75_S_obs")  
# plt.plot(b_t_baseline[2]['time'],state_c_75[0]["S"], label = "75_S")
plt.plot(b_c_15[0]['time'],b_c_15[0]["b(t)"]*3,92*0.05, label = "15")     
plt.xlabel("time") 
plt.ylabel("Transmission rate")
plt.title("Effect of Underreporting and Correction on Suspetible Compartment for policy duration 56 days")
plt.legend()
plt.show()

#%%  
plt.figure(figsize = (10,7)) 
plt.plot(b_t_baseline[0]['time'],b_t_baseline[0]['b(t)'] * 3.92*0.05, label = "duration = 14 days")  
plt.plot(b_t_baseline[1]['time'],b_t_baseline[1]['b(t)'] * 3.92*0.05, label = "duration = 28 days")
plt.plot(b_t_baseline[2]['time'],b_t_baseline[2]['b(t)'] * 3.92*0.05, label = "duration = 56 days")    
plt.xlabel("time") 
plt.ylabel("transmission rate")
plt.legend()
plt.show()
        
        
        
#%%
# calculate the 
#%%
#plot extreme case
file_name_extreme = ['b(t) for case 3 tau = 14 cc = 0.01',
             'b(t) for case 3 tau = 28 cc = 0.01',
             'b(t) for case 3 tau = 56 cc = 0.01',
             ]
b_t_c_extreme = []
for i in range(0,len(file_name_extreme)):
    df = pd.read_csv(outpath + file_name_extreme[i])
    b_t_c_extreme.append(df)

plt.figure(figsize = (10,7)) 
plt.plot(b_t_c_extreme[2]['time'],b_t_c_extreme[2]['b(t)'] * 3.92*0.05, label = "extreme rate 0.01")      
plt.xlabel("time") 
plt.ylabel("transmission rate")
plt.title("Extreme case")
plt.legend()
plt.show()
        
        
        
        
        
        
        
        
        
        
        