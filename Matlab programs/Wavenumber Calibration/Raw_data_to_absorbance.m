%% Raw_data_to_absorbance
% C
% Syntax:
%
%       [
% 
% Input:
%       Reference_Intensity_Raw
%           from 20140703_Data_Processed_Wavenumber_Calibration.mat, which
%           is currently (17-9-2014) located in 'L:\IST\OP\scratch\Olav
%           Grouwstra\' 'Measurements\Wavenumber calibrated data 3july2014'
%       Measurement_Intensity_Raw
%       wavenumber_data
% Output:
%       absorbance_mean_adjusted_data 
%           which is wavenumber calibrated absorbance
%       
%
% Example(s):
%     warning('off','all')
    monitor_intensity_data_ref      =   Processed_Data_Wavenumber_Calibration{1,1}.Reference_Monitor_Intensity_Raw;
    monitor_intensity_data_meas     =   Processed_Data_Wavenumber_Calibration{1,1}.Measurement_Monitor_Intensity_Raw;
    main_intensity_data_ref         =   Processed_Data_Wavenumber_Calibration{1,1}.Reference_MainSignal_Intensity_Raw;
    main_intensity_data_meas        =   Processed_Data_Wavenumber_Calibration{1,1}.Measurement_MainSignal_Intensity_Raw;
    wavenumber_data                 =   Processed_Data_Wavenumber_Calibration{1,1}.Wavenumber;
%     [absorbance_mean_adjusted_data]=Raw_data_to_absorbance(wavenumber_data, monitor_intensity_data_ref, monitor_intensity_data_meas, main_intensity_data_ref, main_intensity_data_meas);
%
% Dependence:   peakfinder
%               Calibration_to_peaks_mean
%
% Author: Olav Grouwstra
% Email: olavgrouwstra@gmail.com
% 
% Project: IR Spectroscopy
% 
%
%Version: v1.0 (21-5-2014)
%
%% Change-log
%
% v1.0 (21-05-2014)  First operation version - Calibration by parts

%% Function declaration
% function [absorbance_mean_adjusted_data]=Raw_data_to_absorbance(wavenumber_data, monitor_intensity_data_ref, monitor_intensity_data_meas, main_intensity_data_ref, main_intensity_data_meas)

close all;

[rowmonref,colmonref]                     =   size(monitor_intensity_data_ref);
[rowmonmeas,colmonmeas]                   =   size(monitor_intensity_data_meas);

[adjusted_intensity_data_monitor,adjusted_intensity_data_main]    =   Calibration_to_peaks_mean(wavenumber_data,[monitor_intensity_data_ref, monitor_intensity_data_meas],[main_intensity_data_ref, main_intensity_data_meas]);

Adjusted_Reference_Intensity            =   adjusted_intensity_data_main(:,1:colmonref)./adjusted_intensity_data_monitor(:,1:colmonref);
Adjusted_Measurement_Intensity          =   adjusted_intensity_data_main(:,colmonref+1:colmonref+colmonmeas)./adjusted_intensity_data_monitor(:,colmonref+1:colmonref+colmonmeas);
Reference_Intensity_Raw                 =   main_intensity_data_ref./monitor_intensity_data_ref;
Measurement_Intensity_Raw               =   main_intensity_data_meas./monitor_intensity_data_meas;
[refrow,refcol]                         =   size(Reference_Intensity_Raw);
[measrow,meascol]                       =   size(Measurement_Intensity_Raw);

Mean_Adjusted_Reference_Intensity               =   mean(Adjusted_Reference_Intensity')';
Mean_Adjusted_Measurement_Intensity             =   mean(Adjusted_Measurement_Intensity')';
Mean_Reference_Intensity_old                    =   mean(Reference_Intensity_Raw')';
Mean_Measurement_Intensity_old                  =   mean(Measurement_Intensity_Raw')';

absorbance_adjusted_data                   =   -log10(Adjusted_Measurement_Intensity./Adjusted_Reference_Intensity);
absorbance_data_old                        =   -log10(Measurement_Intensity_Raw./Reference_Intensity_Raw);

absorbance_mean_adjusted_data                   =   mean(absorbance_adjusted_data,2);
absorbance_mean_data_old                        =   mean(absorbance_data_old,2);


figure;
subplot(2,1,1);plot(wavenumber_data,Adjusted_Reference_Intensity,'r',wavenumber_data,Adjusted_Measurement_Intensity,'b')
legend('Adjusted Reference Intensity Raw','Adjusted Measurement Intensity Raw')
title('Calibrated Refs and Meas')

subplot(2,1,2);plot(wavenumber_data,Mean_Adjusted_Reference_Intensity,'r',wavenumber_data,Mean_Adjusted_Measurement_Intensity,'b')
legend('Reference Intensity','Measurement Intensity')
title('Average of calibrated Refs and Meas')

std_absorbance_data_old     =   std(absorbance_data_old,0,2);
std_absorbance_adjusted_data     =   std(absorbance_adjusted_data,0,2);

figure;
handle_abs_old = plot(wavenumber_data,absorbance_mean_data_old+std_absorbance_data_old,'--b',...
    wavenumber_data,absorbance_mean_data_old-std_absorbance_data_old,'--b',...
    wavenumber_data,absorbance_mean_data_old,'r');
    xlabel('Wavenumber (cm^-^1)')
    ylabel('Absorbance')
legend(handle_abs_old([3 2]),'Mean uncalibrated absorbance','Standard deviation');
title('Mean absorbance of ten measurements')


figure;
handle_abs_new = plot(wavenumber_data,absorbance_mean_adjusted_data+std_absorbance_adjusted_data,'--b',...
    wavenumber_data,absorbance_mean_adjusted_data-std_absorbance_adjusted_data,'--b',...
    wavenumber_data,absorbance_mean_adjusted_data,'r');
    xlabel('Wavenumber (cm^-^1)')
    ylabel('Absorbance')
legend(handle_abs_new([3 2]),'Mean calibrated absorbance','Standard deviation');
title('Mean absorbance of ten measurements')



% Compare ways of getting the std
std_Adjusted_Reference_Intensity                =   std(Adjusted_Reference_Intensity')';
std_Reference_Intensity_Raw                     =   std(Reference_Intensity_Raw')';

std_ratio_Adjusted_Reference_Intensity          =   std((Adjusted_Reference_Intensity(:,1:5)./Adjusted_Reference_Intensity(:,6:10))')';
std_ratio_Reference_Intensity_Raw               =   std((Reference_Intensity_Raw(:,1:5)./Reference_Intensity_Raw(:,6:10))')';

std_log_Adjusted_Reference_Intensity            =   std(-log10(Adjusted_Reference_Intensity)')';
std_log_Reference_Intensity_Raw                 =   std(-log10(Reference_Intensity_Raw)')';

std_log_ratio_Adjusted_Reference_Intensity      =   std(-log10(Adjusted_Reference_Intensity(:,1:5)./Adjusted_Reference_Intensity(:,6:10))')';
std_log_ratio_Reference_Intensity_Raw           =   std(-log10(Reference_Intensity_Raw(:,1:5)./Reference_Intensity_Raw(:,6:10))')';

%   This is a figure showing different representations of the standard
%   deviation of the signal
%     figure;
%     subplot(2,2,1);plot(wavenumber_data,std_Adjusted_Reference_Intensity,'b',wavenumber_data,std_Reference_Intensity_Raw,'r')
%     legend('std Adjusted','std Raw')
%     title('std Reference')
%     subplot(2,2,2);plot(wavenumber_data,std_ratio_Adjusted_Reference_Intensity,'b',wavenumber_data,std_ratio_Reference_Intensity_Raw,'r')
%     legend('std ratio Adjusted','std ratio Raw')
%     title('std ratio Reference')
%     subplot(2,2,3);plot(wavenumber_data,std_log_Adjusted_Reference_Intensity,'b',wavenumber_data,std_log_Reference_Intensity_Raw,'r')
%     legend('std log Adjusted','std log Raw')
%     title('std log Reference')
%     subplot(2,2,4);plot(wavenumber_data,std_log_ratio_Adjusted_Reference_Intensity,'b',wavenumber_data,std_log_ratio_Reference_Intensity_Raw,'r')
%     legend('std log ratio Adjusted','std log ratio Raw')
%     title('std log ratio Reference')


% Prepare std-data
% absorbance_data             =   -log10(Measurement_Intensity_Raw./Reference_Intensity_Raw);
% absorbance_adjusted_data    =   -log10(Adjusted_Measurement_Intensity./Adjusted_Reference_Intensity);
% 
% std_absorbance_data                                 =   std(absorbance_data);
% std_absorbance_adjusted_data                        =   std(absorbance_adjusted_data);
% std_absorbance_calibration_reference_to_intensity   =   std(absorbance_calibration_reference_to_intensity);

% plot();


% end