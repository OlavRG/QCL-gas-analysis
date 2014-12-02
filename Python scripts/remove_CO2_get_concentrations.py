# 1. Get concentrations of molecules by using nonlinear least squares fitting.
# 2. Subtract all absorbances except that of the CO2 molecule from the measurement.
# 3. Fit Voigth peaks to CO2-only spectrum. Now subtract those from total spectrum

import os
import sys
import pandas as pd
import re
import numpy as np
from lmfit import minimize, Parameters

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


# Get absorbance .txt input to dataframes
# Load measurements
breath_spectrum_path = 'Measurements\data of 26-8-2014\data no zeros\\'
breath_spectrum_list = os.listdir(breath_spectrum_path)
for ii in breath_spectrum_list:
    filename = os.path.join(cwd, breath_spectrum_path, ii)
    breath_spectrum = pd.read_csv(filename, sep='\t', header=1).astype(float)


# Get database absorbance .txt input to dataframes
#   interpolation

# Non-linear least squares regression

# Get CO2 peak wavenumbers, intensity, indices

# (Optional peakfinder)

# Fit a gauss to a wavenumber+intensity region

# Again: Non-linear least squares regression
