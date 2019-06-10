

#!/usr/bin/env python
# coding: utf-8


# data analysis and wrangling
import pandas as pd
import numpy as np
import os
from scipy.stats import zscore
from  statistics import stdev



#declare variables
GLM = ("GLM-04")
s = ("01", "02", "03", "04", "05", "06", "07", "09", "10", "11", "12", "13","14", "15", "16", "17","18", "20", "21", "22","23", "24","25", "26")
c = ("control_lik", "neutral_lik", "reward_lik", "control_int", "neutral_int", "reward_int")
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
        if stdev(mod_data[2]) == 0:
            mod_data[2] = 0
            os.chdir(modpath)
            mod_data.to_csv(GLM + '_' + task + '_odor_' + cond + '_fsl.txt', sep='\t', index=False, header=None)
        else:
            mod_data[2] = zscore(mod_data[2])
            mod_data = mod_data.round(5)
            os.chdir(modpath)
            mod_data.to_csv(GLM + '_' + task + '_odor_' + cond + '_fsl.txt', sep='\t', index=False, header=None)
