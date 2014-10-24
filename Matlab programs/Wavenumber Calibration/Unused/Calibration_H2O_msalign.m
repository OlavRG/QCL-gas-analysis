%% Calibration_H2O_msalign
% Input:
%     addpath(genpath('H:/MEP'));
%     load('H:\MEP\Matlab programs\H2O_Acetone_Determination\Database_H2O_and_Acetone.mat')
%     load('H:\MEP\Measurements\MST_Breath_20140328.mat')
%         ref_spectrum_k      =   Database_H2O_and_Acetone.wavenumber;
%         ref_spectrum_int    =   Database_H2O_and_Acetone.H2O_base10_1ppm;
%         meas_data_k         =   MST_Breath_Data{1,1}.Wavenumber;
%         meas_data_int       =   MST_Breath_Data{1,1}.Absorbance;

% Output:
%   
% Dependence: peakfinder
% 
% Author: Olav Grouwstra
% Email: olavgrouwstra@gmail.com
% 
% Project: IR Spectroscopy

%% Function declaration
%Calibrate wavenumbers of measured data with respect to a reference
%molecule
 close all;
    %Define area of interest of the reference spectrum
    ref_spectrum_k_start=1100;
    ref_spectrum_k_end=1300;
    ref_indice_k_start=find(MST_Breath_Data{1,1}.Wavenumber==ref_spectrum_k_start);
    ref_indice_k_end=find(MST_Breath_Data{1,1}.Wavenumber==ref_spectrum_k_end);
    Absorbance_raw= -log10(MST_Breath_Data{1,1}.Measurement_Intensity_Raw...
                    ./MST_Breath_Data{1,1}.Reference_Intensity_Raw);
    Absorbance_raw(isnan(Absorbance_raw))=0;
    Absorbance_raw_ave=mean(Absorbance_raw')';

    %Find peaks to calibrate the measured data to.
    ref_spectrum_int=Absorbance_raw(:,1);
    ref_spectrum_k=MST_Breath_Data{1,1}.Wavenumber;
    [peak_indice,peak_intensity]=peakfinder(ref_spectrum_int);
    peak_k=ref_spectrum_k(peak_indice);
    
    msheatmap(MST_Breath_Data{1,1}.Wavenumber,[Absorbance_raw Absorbance_raw_ave],'markers',peak_k)
    title('before alignment')

%Absorbance_raw_H2O=Absorbance_raw(1170:1350,:)

    YA = msalign(MST_Breath_Data{1,1}.Wavenumber,Absorbance_raw,peak_k,'Rescaling',false);
    msheatmap(MST_Breath_Data{1,1}.Wavenumber,YA,'markers',peak_k)
    title('after alignment')



