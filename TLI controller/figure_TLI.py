#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 14:17:59 2022

@author: gracejia
"""


#%%
# read data 
outpath = '/Users/gracejia/Documents/A-UW//covid19 NSF Project/TLI controller/out/'
import sys
sys.path.append('/Users/gracejia/Documents/A-UW//covid19 NSF Project/TLI controller')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#%%
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
plt.figure(figsize = (16,9)) 
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.plot(b_t_baseline[2]['time'],b_t_baseline[2]['b(t)'] * 3.92*0.05,':', label = "Baseline", color = 'black', lw = 3)  
plt.plot(b_t_nc_1[2]['time'],b_t_nc_1[2]['b(t)'] * 3.92*0.05,'--', label = "Scenario 1 (c = 1)", color = 'green')      
plt.plot(b_t_c_75[2]['time'],b_t_c_75[2]['b(t)'] * 3.92*0.05, label = "Scenario 2 (c = 1.3)", color = 'red',lw = 1)
plt.plot(b_t_c_50[2]['time'],b_t_c_50[2]['b(t)'] * 3.92*0.05, label = "Scenario 3 (c = 2)", color = 'orange')
plt.plot(b_t_c_25[2]['time'],b_t_c_25[2]['b(t)'] * 3.92*0.05, label = "Scenario 4 (c = 4)", color = 'pink')  
plt.plot(b_t_c_10[2]['time'],b_t_c_10[2]['b(t)'] * 3.92*0.05, label = "Scenario 5 (c = 10)", color = 'blue')

plt.xlabel("Days", fontsize = 16) 
plt.ylabel("Transmission Rate", fontsize = 16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# plt.title("Effect of Underreporting and Correction on Transmission rate for policy duration 56 days")
plt.legend(loc = "best", fontsize = 12)
plt.show()

#%%
#plot the number of infectious people (multiply by 100,000)
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

#%%

timeseries = np.linspace(0,360,10000)
plt.figure(figsize = (16,9))
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.plot(timeseries,state_baseline[0]['I'] * 100000,':', label = "Baseline", color = 'black', lw = 3)  
plt.plot(timeseries,state_nc[0]['I']* 100000,'--', label = "Scenario 1 (c = 1)", color = 'green')
plt.plot(timeseries,state_c_75[0]['I']* 100000, label = "Scenario 2 (c = 1.3)", color = 'red',lw = 1)
plt.plot(timeseries,state_c_50[0]['I']* 100000, label = "Scenario 3 (c = 2)", color = 'orange')
plt.plot(timeseries,state_c_25[0]['I']* 100000, label = "Scenario 4 (c = 4)", color = 'pink') 
plt.plot(timeseries,state_c_10[0]['I']* 100000, label = "Scenario 5 (c = 10)", color = 'blue')

plt.xlabel("Days", fontsize = 16) 
plt.ylabel("Number of Infectious People",fontsize = 16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# plt.title("Effect of Underreporting and Correction on Infectious for policy duration 56 days")
plt.legend(loc = "best", fontsize = 12)
plt.show()
# %%
