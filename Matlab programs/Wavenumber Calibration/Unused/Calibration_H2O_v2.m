%% Calibration_H2O_v2
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
% v1.0 (18-02-2014)  First operation version - Global Calibration
% v2.0 (18-02-2014)  Calibration by parts

%% Function declaration
function [ output1 ] = Calibration_H2O_v2(input1, input2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

data_wavenumber = input1;
data_absorbance = input2;

%% Load calibration file
if exist('WavenumberCalibration_H20.mat','file') == 2
    %load the calibration file
    cal_H2O = load('WavenumberCalibration_H20.mat');
    cal_indexes = cal_H2O.H2O_Calibration_peak.indexes;
    cal_wavenumber = cal_H2O.H2O_Calibration_peak.wavenumber;
else
    cal_indexes = 0;
    cal_wavenumber = 0;
end

%% Copy the relevant wavelength region of the data
wavenumber_1150 = find(data_wavenumber == 1150, 1,'first');

data_wavenumber_region = data_wavenumber(wavenumber_1150:end) ;
data_absorbance_region = data_absorbance (wavenumber_1150:end);

size(data_wavenumber_region)
%% Find the closest peak in the data to the calibration peaks
index = [];
for i = 1:length(cal_indexes)
   
    index(i) =  find(data_wavenumber_region >= cal_wavenumber(i,1),1,'first');
    start_search_index = index(i) - 5;
    end_search_index = index(i) + 20;
    search_zone = data_absorbance_region(start_search_index:end_search_index);
    [max_search, max_index]= max(search_zone);
    max_indexes(i) = max_index + start_search_index - 1;
end

data_wavenumber_peaks = data_wavenumber_region(max_indexes);

% figure
% plot(data_wavenumber_region,data_absorbance_region,'r')
% hold on
%% Linear fit between each adjacent peaks.
calibrated_wavenumber = [];
calibrated_absorbance = [];

for j = 1:length(cal_indexes)-1

cuadfit = fit(data_wavenumber_peaks(j:j+1),cal_wavenumber(j:j+1),'poly1');

if j == 1 
    cal_data_wavenumber = feval(cuadfit,data_wavenumber_region(1:max_indexes(j+1)-1));
elseif j == length(cal_indexes)-1
    cal_data_wavenumber = feval(cuadfit,data_wavenumber_region(max_indexes(j):end));
else
    cal_data_wavenumber = feval(cuadfit,data_wavenumber_region(max_indexes(j):max_indexes(j+1)-1));
end

calibrated_wavenumber = [calibrated_wavenumber ; cal_data_wavenumber ];
% calibrated_absorbance = [calibrated_absorbance ; data_absorbance_region(max_indexes(j):max_indexes(j+1)-1) ];

end
 
% size(calibrated_wavenumber)
% size(calibrated_absorbance)

% plot(cal_data_wavenumber,data_absorbance_region(1:max_indexes(j+1)-1))
%% Four order fit

% fitting_ref= cal_wavenumber;
% fitting_meas= data_wavenumber_region(max_indexes);
% cuadfit = fit(fitting_meas,fitting_ref,'poly4');

%  figure;
%  plot (fitting_meas,fitting_ref,'r*');
%  hold on
%  plot (fitting_meas,feval(cuadfit,fitting_meas))

% calibrated_wavenumber = feval(cuadfit,data_wavenumber_region);
% 
 calibrated_data(:,1) = calibrated_wavenumber;
 calibrated_data(:,2) = data_absorbance_region;
% 
 output1 = calibrated_data;

end

