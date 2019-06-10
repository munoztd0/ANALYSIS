#!/usr/bin/env python
# coding: utf-8

"""
Created on Mon Jun 10 14:13:20 2019

@author: david munoz
"""
# data analysis and wrangling
import pandas as pd
import numpy as np
import os


#declare variables
GLM = ("GLM-04")
s = ("01", "02", "03", "04", "05", "06", "07", "09", "10", "11", "12", "13","14", "15", "16", "17","18", "20", "21", "22","23", "24","25", "26")
cond = ("effort") 
taskDIR = ("PIT")

df1 = []
df2 = []
df3 = []
df4 = []
dfsubj = []
dfcontrol = []
df01 = pd.DataFrame()
df02 = pd.DataFrame()
df03 = pd.DataFrame()
df04 = pd.DataFrame()


for i in s:
    subj = 'sub-' + i
    covpath = '/home/cisa/CISA/REWOD/DATA/STUDY/MODELS/SPM/' + taskDIR + '/' + GLM + '/' + subj + '/timing/'
    cov_Base = pd.read_table(covpath + GLM + '_task-PIT_CS_Baseline_fsl.txt',sep='\t', header=None)
    cov_minus = pd.read_table(covpath + GLM + '_task-PIT_CS_CSm_fsl.txt',sep='\t', header=None)
    cov_plus = pd.read_table(covpath + GLM + '_task-PIT_CS_CSp_fsl.txt',sep='\t', header=None)

    dfsubj = np.append(dfsubj, i)

    CSp_CSm = cov_plus[2] - cov_minus[2]
    df1 = np.append(df1, CSp_CSm.mean())


    CSp_Baseline = cov_plus[2] - cov_Base[2]
    df2 = np.append(df2, CSp_Baseline.mean())
    
    
    CSp_CSmandBaseline = cov_plus[2] - (cov_minus[2] + cov_Base[2])/2
    df3 = np.append(df3, CSp_CSmandBaseline.mean())


df01[0] = dfsubj
df02[0] = dfsubj
df03[0] = dfsubj
df01[1] = df1
df02[1] = df2
df03[1] = df3

df01.columns = ['subj', cond]
df02.columns = ['subj', cond]
df03.columns = ['subj', cond]


os.chdir('/home/cisa/CISA/REWOD/DATA/STUDY/MODELS/SPM/PIT/GLM-04/group_covariates')
df01.to_csv('CSp-CSm_' + cond + '.txt',sep='\t', index=False)
df02.to_csv('CSp-Baseline_' + cond + '.txt',sep='\t', index=False)
df03.to_csv('CSp-CSm&Baseline_' + cond + '.txt',sep='\t', index=False)
