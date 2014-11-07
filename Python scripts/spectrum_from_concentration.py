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

if os.path.dirname(os.path.abspath('__file__')) == os.getcwd()
dir = os.path.dirname(os.path.dirname(os.path.abspath('__file__')))
os.chdir(dir)

interaction_length = 54.36
compound_path   =   'L:\IST\OP\scratch\Adonis\Databases\PNNL Database\Compounds\\'
file_regexp     =   '.*_25T?.TXT'
file_extension  =   '.TXT'

file = 'Measurements\data of 26-8-2014\data no zeros\G0S1.txt'
filename = os.path.join(dir, file)
df = pd.read_csv(filename, sep='\t', header=1).astype(float)

