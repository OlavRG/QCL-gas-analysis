# Concentration analysis
# Input
    #G0S1.txt-G2S2.txt
    #wavenumber.txt
    #Compounds from database
    #Directory of standard compounds
    
    
    
# To do:
#- Perform (nonlinear) least squares regression on all compounds (too slow?)
#    - Filter out compounds:
#        - Dump all with concentration<whatever based on single molecule fit

import os
import sys
import pandas as pd
import re
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

# Load standard compounds
wavelength = pd.read_csv(compound_file_path, sep='\t', header=0).astype(float)
interaction_length = 54.36
compound_all_folder   =   '\\Compounds\\'
file_regexp     =   re.compile('.*_25T?.TXT').search
compound_all_list = os.listdir(cwd+compound_all_folder)
for ii in compound_all_list:
    compound_folder = cwd+compound_all_folder+ii
    compound_file_list = os.listdir(compound_folder)
    for l in compound_file_list:
        if file_regexp(l):
            compound_file = file_regexp(l).group(0)
            compound_file_path= os.path.join(compound_folder, compound_file)
            G0S1 = pd.read_csv(compound_file_path, sep='\t', header=0).astype(float)
        #for m in (file_regexp(l),):

        #    print m
            #if m:
            #    print m.group(0)   
 #   compound_path = 5
  #  temporary_compound = pd.read_csv(filename, \
#                                    sep='\t', header=1).astype(float)

# Load measurements
file = 'Measurements\data of 26-8-2014\data no zeros\G0S1.txt'
filename = os.path.join(cwd, file)
G0S1 = pd.read_csv(filename, sep='\t', header=1).astype(float)





def residual(params, x, data, eps_data):
    amp = params['amp'].value
    pshift = params['phase'].value
    freq = params['frequency'].value
    decay = params['decay'].value

    model = amp * sin(x * freq  + pshift) * exp(-x*x*decay)

    return (data-model)/eps_data

params = Parameters()
params.add('amp', value=10)
params.add('decay', value=0.007)
params.add('phase', value=0.2)
params.add('frequency', value=3.0)

out = minimize(residual, params, args=(x, data, eps_data))