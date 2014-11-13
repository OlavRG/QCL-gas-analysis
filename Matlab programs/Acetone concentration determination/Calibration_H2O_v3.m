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
% v2.1 (07-03-2014)  Interpolation of Calibrated absorbance - Adonis Reyes

 
%% Function declaration
function [ output1, output2 ] = Calibration_H2O_v3(input1, input2)
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
wavenumber_1100 = find(data_wavenumber == 1100, 1,'first');


data_wavenumber_region = data_wavenumber(wavenumber_1100:end);
data_absorbance_region = data_absorbance(wavenumber_1100:end);


%% Find the peaks in the data that corresponds to the calibration peaks

length_indexes = length(cal_wavenumber);


for i = 1:length_indexes
    search_index(1,i) =  find(data_wavenumber_region == cal_wavenumber(i,1),1,'first');
    start_search_index = search_index(1,i) - 5;
    end_search_index = search_index(1,i) + 20;
    search_zone = data_absorbance_region(start_search_index:end_search_index);
    
    [t1,t2]=peakfinder(search_zone);
    [max_search1, max_index1] =  max(t2);
    max_index = t1(max_index1);
    max_indexes(i) = max_index + start_search_index - 1;
    
end

data_wavenumber_peaks = data_wavenumber_region(max_indexes);
%plot(data_wavenumber_region(max_indexes),data_absorbance_region(max_indexes),'*k');

% for k = 1:length(max_indexes)-1
%     ind = [];
%     ind = [max_indexes(k)-50 : max_indexes(k)+50];
%     x = data_wavenumber_region(ind);
%     y = data_absorbance_region(ind);
%     f = fit(x,y,'gauss2');
%     [max_fit(k), max_fit_index(k)] = max(feval(f,x));
%     
%     
%     %plot(x(max_fit_index(k)),max_fit(k),'*r');
% 
%     %hold on
%    
%     %plot(f,x,y);
% 
% end
% ratio = max_fit'./cal_intensity(1:end-1);


%% Linear fit between each adjacent peaks.
calibrated_wavenumber = [];
% calibrated_absorbance = [];
calibrated_wavenumber_peaks = [];

for j = 1:length_indexes-1
    
    cuadfit = fit(data_wavenumber_peaks(j:j+1),cal_wavenumber(j:j+1),'poly1');
    cal_data_wavenumber = feval(cuadfit,data_wavenumber_region(max_indexes(j):max_indexes(j+1)-1));
    calibrated_wavenumber = [calibrated_wavenumber ; cal_data_wavenumber ];
    cal_data_wavenumber_peaks = feval(cuadfit,data_wavenumber_region(max_indexes(j)));
    calibrated_wavenumber_peaks = [calibrated_wavenumber_peaks ; cal_data_wavenumber_peaks ];

end

%%Interpolation
wavenumber_interp = [1150 : 0.05 : 1255];%1260.35];
absorbance_cal = data_absorbance_region(max_indexes(1):max_indexes(j+1)-1);
absorbance_interp = interp1(calibrated_wavenumber, absorbance_cal, wavenumber_interp);



%% Output
%calibrated_data(:,1) = wavenumber_interp;
calibrated_data = absorbance_interp;

output1 = wavenumber_interp;
output2 = calibrated_data;

end

