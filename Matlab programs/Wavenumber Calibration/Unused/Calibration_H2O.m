%% Calibration_H2O
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
% v1.0 (18-02-2014)
% First operation version
%

%% Function declaration
function [ output1 ] = Calibration_H2O(input1, input2)
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

%% Copy the relevant region of the database
wavenumber_1150 = find(data_wavenumber == 1150, 1,'first');

data_wavenumber_region = data_wavenumber(wavenumber_1150:end) ;
data_absorbance_region = data_absorbance (wavenumber_1150:end);


%% Find the closest peak in the data to the calibration peaks
index = [];
for i = 1:length(cal_indexes)
   
    index(i) =  find(data_wavenumber_region >= cal_wavenumber(i,1),1,'first');
    start_search_index = index(i) - 5;
    end_search_index = index(i) + 20;
    search_zone = data_absorbance_region(start_search_index:end_search_index);
    max_search= max(search_zone);
    max_index(i) = find(data_absorbance_region == max_search,1,'first');
    
end


%% Four order fit

fitting_ref= cal_wavenumber;
fitting_meas= data_wavenumber_region(max_index);
cuadfit = fit(fitting_meas,fitting_ref,'poly4');

%  figure;
%  plot (fitting_meas,fitting_ref,'r*');
%  hold on
%  plot (fitting_meas,feval(cuadfit,fitting_meas))

calibrated_wavenumber = feval(cuadfit,data_wavenumber_region);

calibrated_data(:,1) = calibrated_wavenumber;
calibrated_data(:,2) = data_absorbance_region;

output1 = calibrated_data;

end

