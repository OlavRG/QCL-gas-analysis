''' This file searches through the compounds in compound_all_folder, and
outputs a .txt file with the folder names of all compounds with any data
above minimum_intensity.'''

import os
import re
import pandas as pd

compound_all_folder = 'Y:\\'
compound_all_list = os.listdir(compound_all_folder)

file_regexp     =   re.compile('.*_25T?.TXT').search
for ii in compound_all_list:
    compound_folder = compound_all_folder+ii
    compound_file_list = os.listdir(compound_folder)
    for l in compound_file_list:
        if file_regexp(l):
            compound_file = file_regexp(l).group(0)
            compound_file_path= os.path.join(compound_folder, compound_file)
            compound_temp = pd.read_csv(compound_file_path, sep='\s', header=0).astype(float)
            compound = compound.append(compound_temp, ignore_index=True)
