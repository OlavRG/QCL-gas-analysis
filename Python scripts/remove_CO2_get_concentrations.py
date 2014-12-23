#!/usr/bin/env python

# 1. Get concentrations of molecules by using nonlinear least squares fitting.
    # Find way to match lsqnonlin to match better with MATLAB
        # Setting upper and lower bounds. This might be done using lmfit-py
# 2. Subtract all absorbances except that of the CO2 molecule from the measurement.
# 3. Fit Voigth peaks to CO2-only spectrum. Now subtract those from total spectrum

# TODO:
# - get script working
# - use while loop to get to right directory in main()

import os
import sys
#import pandas as pd
import re
import numpy as np
import breath_analysis_module as bam
import matplotlib.pyplot as plt

def main():
    # Set the cwd assuming cwd is 'Breath Analysis' or one folder deeper.
    sep = os.sep
    MAIN_PATH_NAME = 'Breath Analysis'
    STR_PATTERN = '.*'+sep+MAIN_PATH_NAME+sep
    MAIN_PATTERN = re.compile(STR_PATTERN.encode('string-escape'))
    
    if  re.match(MAIN_PATTERN,os.getcwd()):
        #os.path.basename(os.getcwd()) == 'Python scripts':
        cwd = os.path.dirname(os.path.dirname(os.path.abspath('__file__')))
        os.chdir(cwd)
    elif os.path.basename(os.getcwd()) == 'Breath Analysis':
        cwd = os.path.dirname(os.path.abspath('__file__'))
    else:
        sys.exit(    "CWD needs to be ~\Breath Analysis\ or at most one folder "
                    "deeper")

if __name__ == '__main__':
    main()
cwd=os.getcwd()


# Load measurements (G0S1.txt) and wavenumbers to arrays
breath_spectrum_folder = 'Measurements\\data of 26-8-2014\\data no zeros\\'
healthy = bam.load_measurement(breath_spectrum_folder, 0, 2)
asthma = bam.load_measurement(breath_spectrum_folder, 1, 2)
abs_meas = healthy[:,0]

#wavenumber = list(np.loadtxt(os.path.join('Measurements','Wavenumber.txt')).astype(str))
wavenumber = np.loadtxt(os.path.join('Measurements','Wavenumber.txt'))
wavenumber = wavenumber[:healthy.shape[0]]


# Load list of compound names 
path_Compounds = 'Compounds'
molecule_list_rel = os.listdir(os.path.join(cwd,path_Compounds))
molecule_list = (os.path.join(cwd,path_Compounds,molecule) for molecule in molecule_list_rel)
molecule_list = filter(os.path.isdir, molecule_list)


file_pattern = '25T.TXT'
file_regexp = re.compile(file_pattern.encode('string-escape'))

# Get database absorbance .txt input to dataframes
standard_compound = []
for folder in molecule_list:
    folder_list = os.listdir(folder)
    if not folder_list:
        molecule_list.remove(folder)
    for file_name in folder_list:
        if re.search(file_regexp,file_name):
            file_path = os.path.join(cwd,path_Compounds,folder,file_name)
            compound_temp = bam.load_database_compound(wavenumber,file_path)
            standard_compound.append(compound_temp)
standard_compound = np.concatenate([standard_compound],1).T
molecule_list_rel = [os.path.relpath(mol,os.path.commonprefix(molecule_list)) for mol in molecule_list]

# Get concentrations using non-linear least squares regression
conc_init = np.ones(standard_compound.shape[1])
gas_length = 54.36
popt, pcov = bam.lsqnonlin(abs_meas,standard_compound,conc_init,gas_length)

# Remove non-CO2 from absorbance
standard_absorbance = standard_compound * popt * gas_length
CO2_index = [re.search('Carbon_dioxide',ind) for ind in molecule_list_rel]
#absorbance_CO2less = abs_meas - 


# Get CO2 peak wavenumbers, intensity, indices

# (Optional peakfinder)

# Fit a gauss to a wavenumber+intensity region

# Again: Non-linear least squares regression

## Plot absorbances
#standard_absorbance = standard_compound * popt * gas_length
#abs_plot = plt.plot(wavenumber, abs_meas, wavenumber, standard_absorbance)
#absorbance_legend = molecule_list_rel
#absorbance_legend.insert(0,'healthy')
#plt.legend(abs_plot, absorbance_legend)
#plt.show()

