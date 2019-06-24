#!/usr/bin/env python
# coding: utf-8




# data analysis and wrangling
import pandas as pd
import os
from scipy.stats import zscore
from  statistics import stdev




#declare variables
GLM = ("GLM-04")
s = ("01", "02", "03", "04", "05", "06", "07", "09", "10", "11", "12", "13","14", "15", "16", "17","18", "20", "21", "22","23", "24","25", "26")
#c = ("CS", "grips")
t = ("Baseline", "CSp", "CSm", "PE", "REM")
task = 'task-PIT'
taskDIR = ("PIT")

for i in s:
    subj = 'sub-' + i
    for j in t:
        cond = j
        # save filepath to variable for easier access
        effortpath = '/home/cisa/CISA/REWOD/DATA/STUDY/MODELS/SPM/' + taskDIR + '/' + GLM + '/' + subj + '/timing/'

        # read the data and store data in DataFrame
        effort_data = pd.read_table(effortpath + GLM + '_' + task + '_CS_' + cond + '.txt',sep='\t', header=None)
        if round(stdev(effort_data[2]) == 0, 7):
            effort_data[2] = 0
            os.chdir(effortpath)
            effort_data.to_csv(GLM + '_' + task + '_CS_' + cond + '_zscored.txt', sep='\t', index=False, header=None)
        else:
            effort_data[2] = zscore(effort_data[2])
            effort_data = round(effort_data, 10)
            os.chdir(effortpath)
            effort_data.to_csv(GLM + '_' + task + '_CS_' + cond + '_zscored.txt', sep='\t', index=False, header=None)
