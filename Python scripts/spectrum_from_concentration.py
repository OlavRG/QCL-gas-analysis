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
import pandas as pd
import re
from lmfit import minimize, Parameters

sep = os.sep
MAIN_PATH_NAME = 'Breath Analysis'
STR_PATTERN = '.*'+sep+MAIN_PATH_NAME+sep
MAIN_PATTERN = re.compile(STR_PATTERN.encode('string-escape')) #'.*'+sep+MAIN_PATH_NAME)

if  re.match(MAIN_PATTERN,os.getcwd()):
    #os.path.basename(os.getcwd()) == 'Python scripts':
    dir = os.path.dirname(os.path.dirname(os.path.abspath('__file__')))
    os.chdir(dir)
elif os.path.basename(os.getcwd()) == 'Breath Analysis':
    dir = os.path.dirname(os.path.abspath('__file__'))
else:
    sys.exit("CWD needs to be ~\Breath Analysis\ or at most one folder deeper")

interaction_length = 54.36
compound_path   =   'L:\IST\OP\scratch\Adonis\Databases\PNNL Database\Compounds\\'
file_regexp     =   '.*_25T?.TXT'
file_extension  =   '.TXT'

file = 'Measurements\data of 26-8-2014\data no zeros\G0S1.txt'
filename = os.path.join(dir, file)
df = pd.read_csv(filename, sep='\t', header=1).astype(float)


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