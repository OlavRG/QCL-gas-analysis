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
breath_spectrum_folder = 'Measurements\\data of 26-8-2014\\data no zeros\\'
healthy = bam.load_measurement(breath_spectrum_folder, 0, 2)
asthma = bam.load_measurement(breath_spectrum_folder, 1, 2)

wavenumber = list(np.loadtxt(os.path.join('Measurements','Wavenumber.txt')).astype(str))
wavenumber = wavenumber[:healthy.shape[0]]


# Load list of compound names 

# Get database absorbance .txt input to dataframes
#   interpolation
molecule_list = np.loadtxt(cwd + '\\Output\\compound_filter\\absorptivity_all_compounds_filtered.txt')
molecule_list = molecule_list[:len(wavenumber),:]

# Non-linear least squares regression
concentration_initial = np.ones(molecule_list.shape[1])
interaction_length = 54.36
bla = bam.lsqnonlin(healthy[:,0],molecule_list,concentration_initial,interaction_length)

# Get CO2 peak wavenumbers, intensity, indices

# (Optional peakfinder)

# Fit a gauss to a wavenumber+intensity region

# Again: Non-linear least squares regression
