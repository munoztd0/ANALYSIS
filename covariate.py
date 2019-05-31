#!/usr/bin/env python
# coding: utf-8


# data analysis and wrangling
import pandas as pd
import os

# %%

#declare variables
GLM = ("GLM-04")
s = ("01", "02", "03", "04", "05", "06", "07", "09", "10", "11", "12", "13","14", "15", "16", "17","18", "20", "21", "22","23", "24","25", "26")
c = ("int") #, "int")
taskDIR = ("hedonic")

df1 = pd.DataFrame()
df2 = pd.DataFrame()
df3 = pd.DataFrame()
df4 = pd.DataFrame()

for i in s:
    subj = 'sub-' + i
    cond = c
    covpath = '/home/cisa/CISA/REWOD/DATA/STUDY/MODELS/SPM/' + taskDIR + '/' + GLM + '/' + subj + '/timing/'
    cov_reward = pd.read_table(covpath + GLM + '_task-hedonic_odor_reward_' + cond + '.txt',sep='\t', header=None)
    cov_control = pd.read_table(covpath + GLM + '_task-hedonic_odor_control_' + cond + '.txt',sep='\t', header=None)
    cov_neutral = pd.read_table(covpath + GLM + '_task-hedonic_odor_neutral_' + cond + '.txt',sep='\t', header=None)
    
    reward_neutral = pd.DataFrame(cov_reward[2] - cov_neutral[2])
    reward_neutral.insert(0, column='subj',value=i)
    reward_neutral.columns = ['subj', cond]  
    df1 = df1.append(reward_neutral)

    reward_control = pd.DataFrame(cov_reward[2] - cov_control[2])
    reward_control.insert(0, column='subj',value=i)
    reward_control.columns = ['subj', cond]
    df2 = df2.append(reward_control)
    
    neutral_control = pd.DataFrame(cov_reward[2] - cov_neutral[2])
    neutral_control.insert(0, column='subj',value=i)
    neutral_control.columns = ['subj', cond]
    df3 = df3.append(neutral_control)
    
    odor_noodor = pd.DataFrame((cov_reward[2] + cov_neutral[2])/2 - cov_control[2])
    odor_noodor.insert(0, column='subj',value=i)
    odor_noodor.columns = ['subj', cond]
    df4 = df4.append(odor_noodor,)
    


os.chdir('/home/cisa/CISA/REWOD/DATA/STUDY/MODELS/SPM/hedonic/GLM-04/group_covariates')
df1.to_csv('reward-neutral_' + cond + '.txt',sep='\t', index=False)
df2.to_csv('reward-control_' + cond + '.txt',sep='\t', index=False)
df3.to_csv('neutral-control_' + cond + '.txt',sep='\t', index=False)
df4.to_csv('odor-noodor_' + cond + '.txt',sep='\t', index=False)
        
        
        
        