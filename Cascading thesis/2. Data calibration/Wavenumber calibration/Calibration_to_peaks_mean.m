 %% Calibration_to_peaks_mean
% C
% Syntax:
%
%       [
% 
% Input:
%       wavenumber_data
%       monitor_intensity_data
%       main_intensity_data
%       optional: begin_range, end_range
% Output:
%   
%       
%
% Example(s):
%       [adjusted_intensity_data]=Calibration_to_peaks_mean(MST_Breath_Data{1,1}.Wavenumber,MST_Breath_Data{1,1}.Measurement_Intensity_Raw(:,:));
%      
% Dependence: peakfinder
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
% wavenumber_data=ref.wavenumber;
% monitor_intensity_data=join_monitor;
% main_intensity_data=join_mainsignal;
% 


%% Function declaration
function [adjusted_intensity_data_monitor,adjusted_intensity_data_main]=Calibration_to_peaks_mean(wavenumber_data,monitor_intensity_data,main_intensity_data,begin_range,end_range)
    
    if exist('begin_range','var')==1 && exist('end_range','var')==1
    else
        begin_range                                 =   wavenumber_data(1);
        end_range                                   =   wavenumber_data(length(wavenumber_data));
    end
        
% Determine the mean by which to calibrate the  peaks and valleys
    [row,col]                                   =   size(monitor_intensity_data);
    mean_intensity                              =   mean(monitor_intensity_data')';
    
% Determine indices of the peaks and valleys of the mean intensity within
% the range determined before.
    index_begin_range                           =   find(wavenumber_data == begin_range,1,'first');
    index_end_range                             =   find(wavenumber_data == end_range,1,'first');        
    peak_limit                                  =   0.0001;
    lower_thresh                                =   0;
    upper_thresh                                =   1.5E9;
    [subindex_peaks_mean,intensity_peak_mean]             =   peakfinder(mean_intensity(index_begin_range:index_end_range),peak_limit,lower_thresh,1);
    [subindex_valleys_mean,intensity_valley_mean]         =   peakfinder(mean_intensity(index_begin_range:index_end_range),peak_limit,upper_thresh,-1);
    index_peaks_mean                            =   index_begin_range + subindex_peaks_mean - 1;
    index_valleys_mean                          =   index_begin_range + subindex_valleys_mean - 1;
    
% Remove front and end peak indices such that area to adjust always starts
% and ends in valleys.
    if index_peaks_mean(1) <= index_valleys_mean(1)
        index_peaks_mean                        =   index_peaks_mean(2:length(index_peaks_mean));
    end
    if index_peaks_mean(length(index_peaks_mean)) >= index_valleys_mean(length(index_valleys_mean))
        index_peaks_mean                        =   index_peaks_mean(1:length(index_peaks_mean)-1);
    end
    
    
    intensity_peak_data_monitor                         =   zeros(length(index_valleys_mean)-1,col);
    subindex_peak_data                          =   zeros(length(index_valleys_mean)-1,col);
    index_peak_data                             =   zeros(length(index_valleys_mean)-1,col);
    intensity_valley_data                       =   zeros(length(index_peaks_mean)-1,col);
    subindex_valley_data                        =   zeros(length(index_peaks_mean)-1,col);
    index_valley_data                           =   zeros(length(index_peaks_mean)-1,col);
    adjusted_intensity_data_monitor             =   monitor_intensity_data;
    adjusted_intensity_data_main                =   main_intensity_data;
    intensity_peak_data_main                    =   zeros(length(index_valleys_mean)-1,col);
    intensity_valley_data_main                  =   zeros(length(index_peaks_mean)-1,col);
    
% Find the values and indices of the peaks of the intensity data. Peaks of
% intensity data can be found inbetween neighbouring valleys,
% whose locations are approximated by the valleys of the mean of the
% measurements. The peaks of the data corresponding to peaks of the mean
% are assumed to be the maximum values found in this region between the
% valleys.
    for k=1:length(index_valleys_mean)-1
        [intensity_peak_data_monitor(k,:),subindex_peak_data(k,:)]  =   max(monitor_intensity_data(index_valleys_mean(k)+1:index_valleys_mean(k+1)-1,:),[],1);
        adjusted_intensity_data_monitor(index_peaks_mean(k),:)      =   intensity_peak_data_monitor(k,:);
        index_peak_data(k,:)                                        =   index_valleys_mean(k) + subindex_peak_data(k,:);
        for m=1:col
            intensity_peak_data_main(k,m)                           =   main_intensity_data(index_peak_data(k,m),m);
        end
        adjusted_intensity_data_main(index_peaks_mean(k),:)         =   intensity_peak_data_main(k,:);   
    end

% Shorten valley mean vector such that minima found will be written to
% those minima valleys that are actually in between peaks within the
% determined range.
    index_valleys_mean                           =  index_valleys_mean(2:length(index_valleys_mean)-1);
    
    for k=1:length(index_peaks_mean)-1
        [intensity_valley_data(k,:),subindex_valley_data(k,:)]  	=   min(monitor_intensity_data(index_peaks_mean(k)+1:index_peaks_mean(k+1)-1,:),[],1);
        adjusted_intensity_data_monitor(index_valleys_mean(k),:)    =   intensity_valley_data(k,:);
        index_valley_data(k,:)                                      =   index_peaks_mean(k) + subindex_valley_data(k,:);
        for m=1:col
            intensity_valley_data_main(k,m)                         =   main_intensity_data(index_valley_data(k,m),m);
        end
        adjusted_intensity_data_main(index_valleys_mean(k),:)       =   intensity_valley_data_main(k,:);   
    end
    
% For each dataset, for each peak/valley to the next valley/peak, the data
% points are projected into the peak-valley space of the adjusted dataset
% where they need to fit, so that using these projections the new data
% points at the wavenumbers corresponding with the peak-valley space of the
% adjusted dataset can be interpolated. This results in
% compressed/stretched datapoints between a peak and the next valley.
    coupled_indices_intensity_peakvalley_data{1,col}           =   [];
    sorted_coupled_indices_intensity_peakvalley_data{1,col}    =   [];
    % Create a sorted list of indices of peaks and valleys as the will
    % appear in the adjusted_intensity_data.
    coupled_indices_intensities_peakvalley                          =   [index_peaks_mean,intensity_peak_data_monitor;index_valleys_mean,intensity_valley_data];
    sorted_coupled_indices_intensity_peakvalley_mean                   =   sortrows(coupled_indices_intensities_peakvalley);

    for k=1:col;
        coupled_indices_intensity_peakvalley_data{1,k}          =   [index_peak_data(:,k),intensity_peak_data_monitor(:,k);index_valley_data(:,k),intensity_valley_data(:,k)];
        sorted_coupled_indices_intensity_peakvalley_data{1,k}   =   sortrows(coupled_indices_intensity_peakvalley_data{1,k});
        for m=1:length(sorted_coupled_indices_intensity_peakvalley_data{1,k})-1
            begin_data                                          =   sorted_coupled_indices_intensity_peakvalley_data{1,k}(m,1);
            end_data                                            =   sorted_coupled_indices_intensity_peakvalley_data{1,k}(m+1,1);
            % Some data points are marked as both a peak and a valley (not
            % sure why), causing a NaN in wavenumber_interval. Following
            % "if" statement circumvents this. 
            if begin_data == end_data
                if m+2 >= length(sorted_coupled_indices_intensity_peakvalley_data{1,k})
                    end_data                                    =   end_data+1;
                else
                    end_data                                        =   sorted_coupled_indices_intensity_peakvalley_data{1,k}(m+2,1);                
                end 
            end
            begin_adjust                                                =   sorted_coupled_indices_intensity_peakvalley_mean(m,1);
            end_adjust                                                  =   sorted_coupled_indices_intensity_peakvalley_mean(m+1,1);
            wavenumber_interval                                         =   (wavenumber_data(end_adjust)-wavenumber_data(begin_adjust))/(end_data-begin_data);
            adjusted_intensity_data_monitor(begin_adjust:end_adjust,k)  =   interp1(wavenumber_data(begin_adjust):wavenumber_interval:wavenumber_data(end_adjust),monitor_intensity_data(begin_data:end_data,k),wavenumber_data(begin_adjust:end_adjust),'linear');
            adjusted_intensity_data_main(begin_adjust:end_adjust,k)     =   interp1(wavenumber_data(begin_adjust):wavenumber_interval:wavenumber_data(end_adjust),main_intensity_data(begin_data:end_data,k),wavenumber_data(begin_adjust:end_adjust),'linear');
        end
    end
                                                    
% Plot old and calibrated intensity data
    figure;
    subplot(2,1,1);%plot(wavenumber_data,mean_intensity,'--rs')
    hold on
    plot(wavenumber_data,monitor_intensity_data(:,[1 2]))
    hold off
    xlabel('Wavenumber (cm^-^1)')
    ylabel('Intensity')
    title('Uncalibrated monitor data of two independent measurements')
    
    subplot(2,1,2);%plot(wavenumber_data,mean_intensity,'--rs')
    hold on
    plot(wavenumber_data,adjusted_intensity_data_monitor(:,[1 2]))
    hold off
    xlabel('Wavenumber (cm^-^1)')
    ylabel('Intensity')
    title('Calibrated monitor data of two independent measurements')

    figure;
    subplot(2,1,1);%plot(wavenumber_data,mean_intensity,'--rs')
    hold on
    plot(wavenumber_data,main_intensity_data(:,[1 2]))
    hold off
    xlabel('Wavenumber (cm^-^1)')
    ylabel('Intensity')
    title('Uncalibrated main data of two independent measurements')
    
    subplot(2,1,2);%plot(wavenumber_data,mean_intensity,'--rs')
    hold on
    plot(wavenumber_data,adjusted_intensity_data_main(:,[1 2]))
    hold off
    xlabel('Wavenumber (cm^-^1)')
    ylabel('Intensity')
    title('Uncalibrated main data of two independent measurements')
end