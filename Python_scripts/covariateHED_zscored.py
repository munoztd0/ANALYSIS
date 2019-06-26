#!/usr/bin/env python
# coding: utf-8

"""
Created on Mon Jun 10 14:13:20 2019

@author: David Munoz
"""
# data analysis and wrangling
import pandas as pd
import numpy as np
import os

# be sure that there is a group_covariates folder

#declare variables
GLM = ("GLM-04")
s = ("01", "02", "03", "04", "05", "06", "07", "09", "10", "11", "12", "13","14", "15", "16", "17","18", "20", "21", "22","23", "24","25", "26")
cond = ("lik") #, "int")
taskDIR = ("hedonic")

df1 = []
df2 = []
df3 = []
df4 = []
df5 = []

dfsubj = []

df01 = pd.DataFrame()
df02 = pd.DataFrame()
df03 = pd.DataFrame()
df04 = pd.DataFrame()
df05 = pd.DataFrame()


for i in s:
    subj = 'sub-' + i
    covpath = '/home/cisa/CISA/REWOD/DATA/STUDY/MODELS/SPM/' + taskDIR + '/' + GLM + '/' + subj + '/timing/'
    cov_control = pd.read_table(covpath + GLM + '_task-hedonic_odor_control_' + cond + '_zscored.txt',sep='\t', header=None)
    cov_neutral = pd.read_table(covpath + GLM + '_task-hedonic_odor_neutral_' + cond + '_zscored.txt',sep='\t', header=None)
    cov_reward = pd.read_table(covpath + GLM + '_task-hedonic_odor_reward_' + cond + '_zscored.txt',sep='\t', header=None)

    dfsubj = np.append(dfsubj, i)

    R_C = cov_reward[2] - cov_control[2]
    df1 = np.append(df1, R_C.mean())

    R_N = cov_reward[2] - cov_neutral[2]
    df2 = np.append(df2, R_N.mean())

    Odor_NoOdor= (cov_reward[2] + cov_neutral[2])/2 - cov_control[2]
    df3 = np.append(df3, Odor_NoOdor.mean())

    Odor_presence = (cov_reward[2] + cov_neutral[2])/2
    df4 = np.append(df4, Odor_presence.mean())

    Reward = cov_reward[2]
    df5 = np.append(df5, Reward.mean())

df01[0] = dfsubj
df02[0] = dfsubj
df03[0] = dfsubj
df04[0] = dfsubj
df05[0] = dfsubj

df01[1] = df1
df02[1] = df2
df03[1] = df3
df04[1] = df4
df05[1] = df5

df01.columns = ['subj', cond]
df02.columns = ['subj', cond]
df03.columns = ['subj', cond]
df04.columns = ['subj', cond]
df05.columns = ['subj', cond]


os.chdir('/home/cisa/CISA/REWOD/DATA/STUDY/MODELS/SPM/hedonic/GLM-04/group_covariates')
df01.to_csv('reward-control_' + cond + '_zscored.txt',sep='\t', index=False)
df02.to_csv('reward-neutral_' + cond + '_zscored.txt',sep='\t', index=False)
df03.to_csv('Odor-NoOdor_' + cond + '_zscored.txt',sep='\t', index=False)
df04.to_csv('Odor_presence_' + cond + '_zscored.txt',sep='\t', index=False)
df05.to_csv('reward_' + cond + '_zscored.txt',sep='\t', index=False)
