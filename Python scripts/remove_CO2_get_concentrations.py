#!/usr/bin/env python


# 1. Get concentrations of molecules by using nonlinear least squares fitting.
    # Find way to match lsqnonlin to match better with MATLAB
        # Setting upper and lower bounds. This might be done using lmfit-py
# 2. Subtract all absorbances except that of the CO2 molecule from the measurement.
    
# 3. Fit Voigth peaks to CO2-only spectrum. Now subtract those from total spectrum

# TODO:
# - get script working
# - use while loop to get to right directory in main()
# - adjust the upkeep of molecule_list and molecule_list_rel

import os
import sys
#import pandas as pd
import re
import numpy as np
import breath_analysis_module as bam
import matplotlib.pyplot as plt
import peakdetect as pdet
from scipy.optimize import curve_fit
from scipy import exp, sqrt

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
cwd = os.getcwd()


# Load measurements (G0S1.txt) and wavenumbers to arrays
breath_spectrum_folder = 'Measurements\\data of 26-8-2014\\data no zeros\\'
healthy = bam.load_measurement(breath_spectrum_folder, 0, 2)
asthma = bam.load_measurement(breath_spectrum_folder, 1, 2)
abs_meas = healthy[:, 0]

# wavenumber = list(np.loadtxt(os.path.join('Measurements','Wavenumber.txt')).astype(str))
wavenumber = np.loadtxt(os.path.join('Measurements', 'Wavenumber.txt'))
wavenumber = wavenumber[:healthy.shape[0]]


# Load list of compound names 
path_Compounds = 'Compounds'
molecule_list_rel = os.listdir(os.path.join(cwd, path_Compounds))
molecule_list = (os.path.join(cwd, path_Compounds, molecule) for molecule in molecule_list_rel)
molecule_list = filter(os.path.isdir, molecule_list)


file_pattern = '25T.TXT'
file_regexp = re.compile(file_pattern.encode('string-escape'))

# Get database absorbance .txt input to dataframes
standard_compound = []
empty_folder = []
for mol_path in molecule_list:
    mol_files = os.listdir(mol_path)
    if not mol_files: #If the mol_files contains no entries (files)
        empty_folder.append(mol_path)
    for file_name in mol_files:
        if re.search(file_regexp,file_name):
            file_path = os.path.join(cwd,path_Compounds,mol_path,file_name)
            compound_temp = bam.load_database_compound(wavenumber,file_path)
            standard_compound.append(compound_temp)
standard_compound = np.concatenate([standard_compound],1).T
{molecule_list.remove(empty) for empty in empty_folder}
molecule_list_rel = [os.path.relpath(mol,os.path.commonprefix(molecule_list)) for mol in molecule_list]

# Get concentrations using non-linear least squares regression
conc_init = np.ones(standard_compound.shape[1])
gas_length = 54.36
conc_opt1, pcov = bam.lsqnonlin(abs_meas, standard_compound, conc_init, gas_length)

# Remove non-CO2 from absorbance
standard_absorbance = standard_compound * conc_opt1 * gas_length
standard_absorbance_noCO2 = standard_absorbance
for idx, val in enumerate(molecule_list_rel):
    if 'Carbon_dioxide' in val:
        standard_absorbance_noCO2 = np.delete(standard_absorbance,idx,1)
        CO2_index = idx
        
abs_meas_onlyCO2 = abs_meas - np.sum(standard_absorbance_noCO2,1)

# Gets CO2 peak using peakdetect, then fits gaussians to them and remove them
abs_remove_all = bam.fit_remove_molecule(abs_meas_onlyCO2, 10, standard_compound[:,CO2_index],wavenumber, 10, 0.0002)

# Add all non CO2 molecules back to the measurement
abs_remove_CO2 = abs_remove_all + np.sum(standard_absorbance_noCO2,1)

# Again: Non-linear least squares regression
conc_opt2, pcov = bam.lsqnonlin(abs_remove_CO2, standard_compound, conc_init, gas_length)

# (Re)state  variables for plotting
abs_meas_onlyCO2 = abs_meas - np.sum(standard_absorbance_noCO2,1) # Restate this variable, somehow fit_remove_molecule changes this variable... 
sum_st_abs = np.sum(standard_absorbance_noCO2,1)

plt.close("all")
plotter = plt.plot(wavenumber,sum_st_abs,wavenumber,abs_meas,wavenumber,abs_meas_onlyCO2, wavenumber, standard_absorbance[:,4], wavenumber, abs_remove_CO2)
plt.legend(plotter, ['summed_noCO2' ,'measurement','meas_onlyCO2', 'DB_CO2', 'abs_remove_CO2'])
plt.show

## Plot absorbances
fig2 = plt.figure()
standard_absorbance = standard_compound * conc_opt1 * gas_length
standard_absorbance2 = standard_compound * conc_opt2 * gas_length
abs_sp1 = fig2.add_subplot(211)
abs_sp2 = fig2.add_subplot(212)
abs_plot1 = abs_sp1.plot(wavenumber, abs_meas, wavenumber, standard_absorbance)
abs_plot2 = abs_sp2.plot(wavenumber, abs_meas, wavenumber, standard_absorbance2)
absorbance_legend = molecule_list_rel
absorbance_legend.insert(0,'healthy')
plt.legend(abs_plot1, absorbance_legend)
#plt.legend(abs_plot2, absorbance_legend)
plt.show()
