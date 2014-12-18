# 1. Get concentrations of molecules by using nonlinear least squares fitting.
# 2. Subtract all absorbances except that of the CO2 molecule from the measurement.
# 3. Fit Voigth peaks to CO2-only spectrum. Now subtract those from total spectrum

# TODO:
# - get script working
# - for performance: use numpy array instead of dataframe
# - use while loop to get to right directory in main()

import os
import sys
import pandas as pd
import re
import numpy as np
import breath_analysis_module as bam

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
wavenumber = list(np.loadtxt('Measurements\Wavenumber.txt').astype(str))
breath_spectrum = np.empty((len(wavenumber),0))
breath_spectrum_folder = 'Measurements\data of 26-8-2014\data no zeros\\'
breath_spectrum_file_list = os.listdir(breath_spectrum_folder)
spectrum_shape = ()
for ii in breath_spectrum_file_list:
    breath_spectrum_file_path = os.path.join(cwd,breath_spectrum_folder,ii)
    breath_spectrum_temp = np.loadtxt(breath_spectrum_file_path, skiprows=1) # Later append requires transpose
    breath_spectrum = breath_spectrum[:len(breath_spectrum_temp)] if (ii==breath_spectrum_file_list[0]) else breath_spectrum
    spectrum_shape = spectrum_shape + breath_spectrum_temp.shape
    breath_spectrum = np.append(breath_spectrum, breath_spectrum_temp, axis=1)

healthy_col1 = spectrum_shape[1] + spectrum_shape[3]
healthy_col2 = healthy_col1 + spectrum_shape[5] + spectrum_shape[7]
asthma_col = healthy_col2 + spectrum_shape[9] + spectrum_shape[11]
healthy = breath_spectrum[:,:healthy_col1]
asthma = breath_spectrum[:,healthy_col2:asthma_col]

# Load list of compound names 

# Get database absorbance .txt input to dataframes
#   interpolation
# 
molecule_list = np.loadtxt(cwd + '\\Output\\compound_filter\\absorptivity_all_compounds_filtered.txt')
molecule_list = molecule_list[:len(breath_spectrum_temp),:]

# Non-linear least squares regression
concentration_initial = np.ones(molecule_list.shape[1])
interaction_length = 54.36
bla = bam.lsqnonlin(healthy[:,0],molecule_list,concentration_initial,interaction_length)

# Get CO2 peak wavenumbers, intensity, indices

# (Optional peakfinder)

# Fit a gauss to a wavenumber+intensity region

# Again: Non-linear least squares regression
