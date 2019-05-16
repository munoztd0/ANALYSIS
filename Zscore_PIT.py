#!/usr/bin/env python
# coding: utf-8

# In[8]:



#!/usr/bin/env python
# coding: utf-8


# data analysis and wrangling
import pandas as pd
import numpy as np
from scipy.stats import zscore
import os


# In[33]:



#declare variables
GLM = ("GLM-03")
s = ("01", "02", "03", "04", "05", "06", "07", "09", "10", "11", "12", "13","14", "15", "16", "17","18", "20", "21", "22","23", "24","25", "26") 
c = ("CSp", "CSm", "Baseline")
t = ("task-PartialExtinction_trial", "task-Reminder_trial")
task = 'task-PIT'
taskDIR = ("PIT")

for i in s:
    subj = 'sub-' + i
    for j in c:
        cond = j
        # save filepath to variable for easier access
        effortpath = '/home/cisa/CISA/REWOD/DATA/STUDY/MODELS/SPM/' + taskDIR + '/' + GLM + '/' + subj + '/timing/' 

        # read the data and store data in DataFrame
        effort_data = pd.read_table(effortpath + GLM + '_' + task + '_CS_' + cond + '.txt',sep='\t', header=None)
        effort_data.iloc[:,2] = zscore(effort_data.iloc[:,2])
        os.chdir(effortpath)
        effort_data.to_csv(GLM + '_' + task + '_CS_' + cond + '.txt', sep='\t', index=False, header=None)

#now with reminder and partial ext
for i in s:
    subj = 'sub-' + i
    for k in t:
        taskP = k
        effort_data = pd.read_table(effortpath + GLM + '_' + taskP + '.txt', sep='\t', header=None)
        effort_data.iloc[:,2] = zscore(effort_data.iloc[:,2])   
        effort_data.to_csv(GLM + '_' + taskP + '.txt', sep='\t', index=False, header=None)



