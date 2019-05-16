#!/usr/bin/env python
# coding: utf-8

# In[1]:



#!/usr/bin/env python
# coding: utf-8


# data analysis and wrangling
import pandas as pd
import numpy as np
from scipy.stats import zscore
import os


# In[3]:



#declare variables
GLM = ("GLM-03")
s = ("01", "02") #, "03", "04", "05", "06", "07", "09", "10", "11", "12", "13","14", "15", "16", "17","18", "20", "21", "22","23", "24","25", "26") 
c = ("control", "neutral", "reward", "int_control", "int_neutral", "int_reward")
task = 'task-hedonic'
taskDIR = ("hedonic")

for i in s:
    subj = 'sub-' + i
    for j in c:
        cond = j
        # save filepath to variable for easier access
        modpath = '/home/cisa/CISA/REWOD/DATA/STUDY/MODELS/SPM/' + taskDIR + '/' + GLM + '/' + subj + '/timing/' 

        # read the data and store data in DataFrame
        mod_data = pd.read_table(modpath + GLM + '_' + task + '_odor_' + cond + '.txt',sep='\t', header=None)
        mod_data.iloc[:,2] = zscore(mod_data.iloc[:,2])
        os.chdir(modpath)
        mod_data.to_csv(GLM + '_' + task + '_odor_' + cond + '.txt', sep='\t', index=False, header=None)


# In[ ]:




