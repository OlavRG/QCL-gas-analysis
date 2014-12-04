# 

import numpy as np

# load_database_compound loads database compounds and fits them to wavenumber
def load_database_compound(wavenumber, compound_path):
    compound1 = np.loadtxt(compound_path)
    
    compound_start_index     = np.argmin(abs(compound1[:,0]-wavenumber[0]))
    compound_end_index     = np.argmin(abs(compound1[:,0]-wavenumber[-1]))

    # This 'if' is used to flip the compound data vectors, which
    # tend to have the highest k in (1,1), and the lowest in (end) for PNNL
    # data.
    if compound_start_index > compound_end_index:
        compound2   =   np.flipud(compound1)
        lower_compound_index     =   len(compound2)-compound_start_index-1-1
        upper_compound_index     =   len(compound2)-compound_end_index-1+1
    else:
        compound2   =   compound1
        lower_compound_index     =   compound_start_index-1
        upper_compound_index     =   compound_end_index+1
    

    compound3     =   interp1( 
    compound2[lower_compound_index:upper_compound_index,1], 
    compound2[lower_compound_index:upper_compound_index,2], 
    wavenumber)

    
    return compound_end_index
    
    
if __name__ == '__main__':
    wavenumber = np.loadtxt('D:\Workspace\Breath Analysis\Measurements\Wavenumber.txt').astype(float)
    compound_path = 'D:\Workspace\Breath Analysis\Compounds\Acetone\ACETONE_25T.TXT'
    bla = load_database_compound(wavenumber,compound_path)