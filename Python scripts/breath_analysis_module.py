#!/usr/bin/env python


import profile
import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
import peakdetect as pdet

def load_measurement(measurement_folder, group, sample_max):
    #group = 0
    #sample_max = 2
    measurement_folder = "D:\\Workspace\\Breath Analysis\\Measurements\\data of 26-8-2014\\data no zeros\\"
    measurement = []
    for smpl in range(1, sample_max+1):
        measurement_file = measurement_folder + 'G%dS%d.txt' % (group, smpl)
        measurement_temp = np.loadtxt(measurement_file, skiprows=1)
        measurement.append(measurement_temp)
    measurement = np.concatenate(measurement,1)
    
    return measurement

# load_database_compound loads database compounds and fits them to wavenumber
# Gives exact same result as my load_compound.m up to 16th significant figure, 
# however this function sometimes has 16 and sometimes 17 significant figures, 
# whereas the MATLAB version always has 16.
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
    
    
    
def lsqnonlin(absorbance, absorptivity_database_molecule_all, concentration_initial, interaction_length):
    def func(absorptivity_length, *C):
        return np.dot(absorptivity_length,C)
    xdata = absorptivity_database_molecule_all * interaction_length
    ydata = absorbance
    popt, pcov = curve_fit(func, xdata, ydata, concentration_initial)
    return popt, pcov


def fit_remove_molecule(absorbance, peakwidth, mlcl, wavenumber, look, Delta):
    mlcl_peaks, mlcl_valleys = pdet.peakdetect(mlcl, lookahead = look, delta = Delta)
    peak_index, peak_height = zip(*mlcl_peaks)
    
    def gauss_func(peak_x, a, x0, sigma):
        return a*np.exp(-(peak_x-x0)**2/(2*sigma**2))
    
    for ind in peak_index:
        peak_range = range(ind-peakwidth,ind+peakwidth)
        peak_x = wavenumber[peak_range]
        peak_y = absorbance[peak_range]
        if max(peak_y) >= 0.04:
            len_x = len(peak_x)
            mean = wavenumber[ind]
            sigma = np.sqrt(np.dot(peak_x-mean,peak_x-mean)/len_x)
            popt, pcov = curve_fit(gauss_func, peak_x, peak_y, p0 = [1, mean, sigma])
            y_fit = gauss_func(peak_x, popt[0], popt[1], popt[2])
            absorbance[ind-peakwidth:ind+peakwidth] = absorbance[ind-peakwidth:ind+peakwidth] - y_fit
        
        return absorbance

#if __name__ == '__main__':
# Test for load_database_compound
    #wavenumber = np.loadtxt('D:\Workspace\Breath Analysis\Measurements\Wavenumber.txt').astype(float)
    #compound_path = 'D:\Workspace\Breath Analysis\Compounds\Acetone\ACETONE_25T.TXT'
    ##bla = profile.run('load_database_compound(wavenumber,compound_path)') #
    #bla = load_database_compound(wavenumber,compound_path)#

# Test lsqnonlin
# bla2 = lsqnonlin()

#Test load_measurement
    #group = 0
    #sample_max = 2
    #measurement_folder = "D:\\Workspace\\Breath Analysis\\Measurements\\data of 26-8-2014\\data no zeros\\"
    #bla = load_measurement(measurement_folder, group, sample_max)