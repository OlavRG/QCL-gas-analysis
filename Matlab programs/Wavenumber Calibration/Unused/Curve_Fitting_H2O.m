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
function [ output1, output2 ] = Curve_Fitting_H2O(input1, input2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

data_wavenumber = input1;
data_absorbance = input2;

%% Load calibration file

if exist('WavenumberCalibration_H20.mat','file') == 2
    %load the calibration file
    cal_H2O = load('WavenumberCalibration_H20.mat');
    cal_wavenumber = cal_H2O.H2O_Calibration_peaks.wavenumber;
    
    cal_intensity = cal_H2O.H2O_Calibration_peaks.intensity;
else
    cal_wavenumber = 0;
end

%% Copy the relevant wavelength region of the data
%wavenumber_1100 = find(data_wavenumber == 1100, 1,'first')


data_wavenumber_region = data_wavenumber;
data_absorbance_region = data_absorbance;


%% Find the peaks in the data that corresponds to the calibration peaks

length_indexes = length(cal_wavenumber);


for i = 2:length_indexes-1
    search_index(1,i) =  find(data_wavenumber_region == cal_wavenumber(i,1),1,'first');
    start_search_index = search_index(1,i) - 5;
    end_search_index = search_index(1,i) + 10;
    search_zone = data_absorbance_region(start_search_index:end_search_index);
    
    [t1,t2]=peakfinder(search_zone);
    [max_search1, max_index1] =  max(t2);
    max_index = t1(max_index1);
    max_indexes(i) = max_index + start_search_index - 1;
    
end

%data_wavenumber_peaks = data_wavenumber_region(max_indexes);
%plot(data_wavenumber_region(max_indexes),data_absorbance_region(max_indexes),'*k');

for k = 2:length(max_indexes)
    ind = [];
    ind = [max_indexes(k)-50 : max_indexes(k)+50];
    x = data_wavenumber_region(ind);
    y = data_absorbance_region(ind);
    f = fit(x,y,'gauss2');
    curve_fitting{k}(:,1) = x;
    y_fitting = feval(f,x);
    curve_fitting{k}(:,2) = y_fitting;
    [max_fit(k), max_fit_index(k)] = max(y_fitting);
    
    
    %plot(x(max_fit_index(k)),max_fit(k),'*r');

    %hold on
   
    %plot(f,x,y);

end


ratio = max_fit'./cal_intensity(2:end);


%% Output

output1 = curve_fitting;
output2 = ratio;

end

