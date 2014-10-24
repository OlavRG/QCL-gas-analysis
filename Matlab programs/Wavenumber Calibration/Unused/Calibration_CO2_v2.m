%% Calibration_CO2_v2
% C
% Syntax:
%
%       [
% 
% Input:
%         
%       
%               
% Output:
%   
%       
%
% Example(s):
%        
%      
% Dependence: peakfinder
% 
% Author: Adonis Reyes Reyes
% Email: adonisatl@gmail.com
% 
% Project: IR Spectroscopy
% 
%
%Version: v1.0 (18-02-2014)
%
%% Change-log
%
% v1.0 (20-02-2014)  First operation version - Calibration by parts
% v2.0 (05-03-2014)  Calibration using peakfinder - Adonis Reyes
% v2.1 (07-03-2014)  Interpolation of Calibrated absrobance - Adonis Reyes

 
%% Function declaration
function [ output1 ] = Calibration_CO2_v2(input1, input2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

data_wavenumber = input1;
data_absorbance = input2;

%% Load calibration file

if exist('WavenumberCalibration_CO2.mat','file') == 2
    %load the calibration file
    cal_CO2 = load('WavenumberCalibration_CO2.mat');
    cal_indexes = cal_CO2.CO2_Calibration_peaks.indexes;
    cal_wavenumber = cal_CO2.CO2_Calibration_peaks.wavenumber;
    cal_intensity = cal_CO2.CO2_Calibration_peaks.intensity;
else
    cal_indexes = 0;
    cal_wavenumber = 0;
end


%% Copy the relevant wavelength region of the data
wavenumber_900 = find(data_wavenumber == 900, 1,'first');
wavenumber_1150 = find(data_wavenumber == 1150, 1,'first');

data_wavenumber_region = data_wavenumber(wavenumber_900:wavenumber_1150) ;
data_absorbance_region = data_absorbance(wavenumber_900:wavenumber_1150);


%% Find the peaks in the data that corresponds to the calibration peaks

length_indexes = length(cal_indexes);
%search_index = zeros(length_indexes);

for i = 1:length_indexes
    search_index(1,i) =  find(data_wavenumber_region == cal_wavenumber(i,1),1,'first');
    start_search_index = search_index(1,i) - 5;
    end_search_index = search_index(1,i) + 10;
    search_zone = data_absorbance_region(start_search_index:end_search_index);
    
    [t1,t2]=peakfinder(search_zone);
    [max_search1, max_index1] =  max(t2);
    max_index = t1(max_index1);
    max_indexes(i) = max_index + start_search_index - 1;

end

data_wavenumber_peaks = data_wavenumber_region(max_indexes);

%% Linear fit between each adjacent peaks.
calibrated_wavenumber = [];
calibrated_absorbance = [];

for j = 1:length(cal_indexes)-1
    
    cuadfit = fit(data_wavenumber_peaks(j:j+1),cal_wavenumber(j:j+1),'poly1');
    cal_data_wavenumber = feval(cuadfit,data_wavenumber_region(max_indexes(j):max_indexes(j+1)-1));
    calibrated_wavenumber = [calibrated_wavenumber ; cal_data_wavenumber ];

end
 
%%Interpolation
wavenumber_interp = [900 : 0.05 : 1150];
absorbance_cal = data_absorbance_region(max_indexes(1):max_indexes(j+1)-1);
absorbance_interp = interp1(calibrated_wavenumber, absorbance_cal, wavenumber_interp);

%% Output
calibrated_data(:,1) = wavenumber_interp;
calibrated_data(:,2) = absorbance_interp;

output1 = calibrated_data;

end

