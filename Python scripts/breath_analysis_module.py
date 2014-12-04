# 

import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit


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
        lower_compound_index     =   len(compound2)-compound_start_index-1-1     # -1 because len is no index
        upper_compound_index     =   len(compound2)-compound_end_index-1+1 
    else:
        compound2   =   compound1
        lower_compound_index     =   compound_start_index-1
        upper_compound_index     =   compound_end_index+1
    

    compound3     =   interpolate.interp1d(
    compound2[lower_compound_index:upper_compound_index+1,0], 
    compound2[lower_compound_index:upper_compound_index+1,1])(wavenumber)

    
    return compound3
    
    
    
def lsqnonlin(absorbance, molecule_list, concentration_initial, interaction_length):
    def func(absorptivity_length, C):
        return absorptivity_length * C
    
    load_database_compound(wavenumber, compound_path)
    xdata = absorptivity_database_compound_all * interaction_length
    ydata = absorbance
    popt, pcov = curve_fit(func, xdata, ydata, concentration_initial)

if __name__ == '__main__':
    #wavenumber = np.loadtxt('D:\Workspace\Breath Analysis\Measurements\Wavenumber.txt').astype(float)
    #compound_path = 'D:\Workspace\Breath Analysis\Compounds\Acetone\ACETONE_25T.TXT'
    #bla = load_database_compound(wavenumber,compound_path)
    
    molecule_list =   list(np.loadtxt('Measurements\Wavenumber.txt').astype(str))