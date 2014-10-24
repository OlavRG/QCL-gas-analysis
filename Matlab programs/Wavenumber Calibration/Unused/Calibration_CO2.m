%% Calibration_CO2
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
% Dependence: NONE (no additional special function are required)
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
 
%% Function declaration
function [ output1 ] = Calibration_CO2(input1, input2)
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
else
    cal_indexes = 0;
    cal_wavenumber = 0;
end

%% Copy the relevant wavelength region of the data
wavenumber_900 = find(data_wavenumber == 900, 1,'first');
wavenumber_1150 = find(data_wavenumber == 1150, 1,'first');

data_wavenumber_region = data_wavenumber(wavenumber_900:wavenumber_1150) ;
data_absorbance_region = data_absorbance(wavenumber_900:wavenumber_1150);

% size(data_wavenumber_region)
% size(data_absorbance_region)
%% Find the closest peak in the data to the calibration peaks
index = [];
for i = 1:length(cal_indexes)
   index(i) =  find(data_wavenumber_region >= cal_wavenumber(i,1),1,'first');
    start_search_index = index(i) - 5;
    end_search_index = index(i) + 15;
    search_zone = data_absorbance_region(start_search_index:end_search_index);
    [max_search, max_index]= max(search_zone);
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
 
% plot(cal_data_wavenumber,data_absorbance_region(1:max_indexes(j+1)-1))
%% Four order fit

calibrated_data(:,1) = calibrated_wavenumber;
calibrated_data(:,2) = data_absorbance_region(max_indexes(1):max_indexes(j+1)-1);
% 
 output1 = calibrated_data;

end

