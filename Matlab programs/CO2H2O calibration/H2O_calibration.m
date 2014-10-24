% Input:
% 
% Output:
%   
%
% Example(s):
%      
% Dependence: peakfinder
% 
% Author: Olav Grouwstra
% Email: olavgrouwstra@gmail.com
% 
% Project: IR Spectroscopy
%
%
% Change-log Version: 
% 
%

% Function declaration
function [adjustedH2O_absorbance,index_peaks_H2O,index_peak_data]=H2O_calibration(wavenumber,absorbance,H2O,index_lastCO2peak)

%     absorbance=adjustedCO2_absorbance;
%     wavenumber  =   Processed_Data_Wavenumber_Calibration{1,1}.Wavenumber;
%     H2O_PNNL  =   flipud(importdata(['L:\IST\OP\scratch\Olav Grouwstra\Compounds\' ...
%         'H2O_25T.TXT']));   
%     H2O=interp1(H2O_PNNL(:,1),H2O_PNNL(:,2),wavenumber);
%     index_lastCO2peak   =   5239;

tic;
    [row,col]                                   =   size(absorbance);

% Determine indices of the peaks of H2O, and identify the max data point
% within index_range of each peak. Now determine the mean distance of the data peaks to the
% H2O peaks, and shift
    peak_limit                                  =   0.01;
    lower_thresh                                =   0.02;
    H2O_factor                                  =   50000; %1ppm concentration of CO2 is too small for calibration
    [index_peaks_H2O,intensity_peak_H2O]        =   peakfinder(H2O*H2O_factor,peak_limit,lower_thresh,1);

    index_peak_data         = zeros(length(index_peaks_H2O),col);
    shiftedH2O_absorbance  = zeros(row,col);
    
    index_range     =   35;
    for k=1:length(index_peaks_H2O)
        [intensity_peak_data,subindex_peak_data]       =   max(absorbance(index_peaks_H2O(k)-index_range:index_peaks_H2O(k)+index_range,:));
        index_peak_data(k,:)    =   subindex_peak_data + index_peaks_H2O(k)-index_range-1;
    end
    mean_peak_error         =   round(mean(mean(index_peak_data,2)-index_peaks_H2O));
    shiftedH2O_absorbance(1:end-mean_peak_error+1,:)  =   [absorbance(1:index_lastCO2peak+1-mean_peak_error,:); ...
        absorbance(index_lastCO2peak+1:end,:)]; 
    shifted_index_peak_data         =   index_peak_data-mean_peak_error+1;
    

% 
        adjustedH2O_absorbance  = shiftedH2O_absorbance;
    for j=1:col
        for k=1:length(index_peaks_H2O)
            begin_data1             =   index_peaks_H2O(k)-index_range-mean_peak_error;
            end_data1               =   shifted_index_peak_data(k,j);
            begin_adjust1           =   index_peaks_H2O(k)-index_range-mean_peak_error;
            end_adjust1             =   index_peaks_H2O(k);
            wavenumber_interval1    =   (wavenumber(end_adjust1)-wavenumber(begin_adjust1))/(end_data1-begin_data1);
            adjustedH2O_absorbance(begin_adjust1:end_adjust1,j)  =   interp1(wavenumber(begin_adjust1):wavenumber_interval1:wavenumber(end_adjust1),shiftedH2O_absorbance(begin_data1:end_data1,j),wavenumber(begin_adjust1:end_adjust1),'linear');

            begin_data2             =   shifted_index_peak_data(k,j);
            end_data2               =   index_peaks_H2O(k)+index_range;
            begin_adjust2           =   index_peaks_H2O(k);
            end_adjust2             =   index_peaks_H2O(k)+index_range;
            wavenumber_interval2    =   (wavenumber(end_adjust2)-wavenumber(begin_adjust2))/(end_data2-begin_data2);
            adjustedH2O_absorbance(begin_adjust2:end_adjust2,j)  =   interp1(wavenumber(begin_adjust2):wavenumber_interval2:wavenumber(end_adjust2),shiftedH2O_absorbance(begin_data2:end_data2,j),wavenumber(begin_adjust2:end_adjust2),'linear');

        end
    end

% absplot     =   2;
% close all
% figure
% hold on
% plot(wavenumber,shiftedH2O_absorbance(:,absplot), ...
%     wavenumber, adjustedH2O_absorbance(:,absplot))
% plot(wavenumber, H2O*H2O_factor,'r')
% plot(wavenumber(index_peaks_H2O),H2O(index_peaks_H2O)*H2O_factor,'xk')
% plot(wavenumber(shifted_index_peak_data(:,absplot)),shiftedH2O_absorbance(shifted_index_peak_data(:,absplot),absplot),'xk')
% hold off
% legend('adjusted absorbance','adjustedH2O_absorbance','H2O','peaks H2O','peaks data')

% toc
end