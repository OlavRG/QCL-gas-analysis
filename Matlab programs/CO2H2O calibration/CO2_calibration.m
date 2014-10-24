% Input:
%       G0S1 etc. matrices of absorbances, w/o sample names (e.g. 022/03)
%       H2O spectrum of PNNL
%       CO2 spectrum of HITRAN
% Output:
%       adjustedH2OCO2_absorbance 
%
% Example(s):
%      
% Dependence: peakfinder
%             H2O_calibration
%             WriteCalibratedDataToText
% 
% Author: Olav Grouwstra
% Email: olavgrouwstra@gmail.com
% 
% Project: IR Spectroscopy
%
% Version: v1.1 (22-5-2014)
%
% Change-log Version: 
% 
% v1.1 (22-5-2014): Now starts looking for next peak/valley from last
% valley/peak found. Works for all of the CO2 peaks 
%
% v1.0 (21-5-2014): Look for data peaks/valleys between two CO2
% valleys/peaks. Works for ~90% of the peaks.
%


% Function declaration
%function [adjusted_absorbance]=CO2_calibration(wavenumber,absorbance,CO2)
tic;

sampleControlPath=['L:\IST\OP\scratch\Olav Grouwstra\Measurements\'...
    'CO2H2O calibrated data\data of 26-8-2014\'...
    'Sample Control_Checked_20140826_with Removed entries.xlsx'];

saveTXTsTo=['L:\IST\OP\scratch\Olav Grouwstra\Measurements\'...
            'CO2H2O calibrated data\data of 26-8-2014\data\'];

% Create G0S1 and so forth using SortSamples.m
    absorbance  =   [G0S1,G0S2,G0S3,G0S4,G1S1,G1S2,G2S1,G2S2];
    [row,col]                                   =   size(absorbance);
    wavenumber=importdata('L:\IST\OP\scratch\Olav Grouwstra\Measurements\Wavenumber.txt');
    H2O_PNNL  =   flipud(importdata(['L:\IST\OP\scratch\Olav Grouwstra\Compounds\' ...
        'H2O_25T.TXT']));   
    H2O=interp1(H2O_PNNL(:,1),H2O_PNNL(:,2),wavenumber);
    CO2_HITRAN  =   importdata(['L:\IST\OP\scratch\Olav Grouwstra\Compounds\' ...
        'CO2_HITRAN.TXT']);   
    CO2 = interp1(CO2_HITRAN.data(:,1), CO2_HITRAN.data(:,2), wavenumber);

toc
    
% Determine indices of the peaks and valleys of CO2
    peak_limit                                  =   0.01;
    lower_thresh                                =   0;
    upper_thresh                                =   1.5E9;
    CO2_factor                                  =   60; %1ppm concentration of CO2 is too small for calibration
    [index_peaks_CO2,intensity_peak_CO2]        =   peakfinder(CO2*CO2_factor,peak_limit,lower_thresh,1);
    [index_valleys_CO2,intensity_valley_CO2]    =   peakfinder(CO2*CO2_factor,peak_limit,upper_thresh,-1);

% Remove front and end peak indices such that area to adjust always starts
% and ends in valleys.
    if index_peaks_CO2(1) <= index_valleys_CO2(1)
    index_peaks_CO2                        =   index_peaks_CO2(2:length(index_peaks_CO2));
    end
    if index_peaks_CO2(length(index_peaks_CO2)) >= index_valleys_CO2(length(index_valleys_CO2))
        index_peaks_CO2                        =   index_peaks_CO2(1:length(index_peaks_CO2)-1);
    end

    intensity_peak_data_absorbance              =   zeros(length(index_valleys_CO2)-1,col);
    subindex_peak_data                          =   zeros(length(index_valleys_CO2)-1,col);
    index_peak_data                             =   zeros(length(index_valleys_CO2)-1,col);
    intensity_valley_data                       =   zeros(length(index_peaks_CO2)-1,col);
    subindex_valley_data                        =   zeros(length(index_peaks_CO2)-1,col);
    index_valley_data                           =   zeros(length(index_peaks_CO2)-1,col);
    adjustedCO2_absorbance                      =   absorbance;
    col_length                                  =   zeros(col,1);
    index_last_valley                           =   zeros(1,col);
    absorbance_matrix                           =   zeros(row,col);

% Look for the first peak between the first two valleys of CO2. Then find
% the next valleys/peaks by looking between the last peak/valley found and
% the next peak/valley from CO2, alternating between finding a peak and a
% valley. This way all valleys are found between two peaks, and all peaks
% are found between two valleys. This loop provides the positions of the
% absorbance peaks/valleys in index_peak_data/index_valley_data.

    % Shorten valley CO2 vector such that minima found will be written to
    % those minima valleys that are actually in between peaks within the
    % determined range.
        index_valleys_CO2_holder                           =  index_valleys_CO2(2:length(index_valleys_CO2)-1);
    for k=1:length(index_valleys_CO2)-1
        if k==1
            index_last_valley(1:col) = index_valleys_CO2(1)+1;
        else
            index_last_valley = index_valley_data(k-1,:);
        end
        
        absorbance_matrix(1:row,1:col)              =   NaN;
        col_length(1:col)  =   index_valleys_CO2(k+1)-1 - index_last_valley+1;
        for j=1:col
            absorbance_matrix(1:col_length(j),j)  = ...
                absorbance(index_last_valley(j):index_valleys_CO2(k+1)-1,j);
        end
        [intensity_peak_data_absorbance(k,:),subindex_peak_data(k,:)]   =   max(absorbance_matrix,[],1);
        index_peak_data(k,:)                                            =   index_last_valley + subindex_peak_data(k,:) - 1;


        if k==length(index_peaks_CO2) %Last part of loop (the valley determination) goes untill length(index_peaks_CO2)=(length(index_valleys_CO2)-2)
            break
        end
        
        index_last_peak  =   index_peak_data(k,:);
        absorbance_matrix(1:row,1:col)              =   NaN;
        col_length(1:col)  =   index_peaks_CO2(k+1)-1 - index_last_peak+1;
        for j=1:col
            absorbance_matrix(1:col_length(j),j)  = ...
                absorbance(index_last_peak(j):index_peaks_CO2(k+1)-1,j);
        end

        [intensity_valley_data(k,:),subindex_valley_data(k,:)]  	=   min(absorbance_matrix,[],1);
        index_valley_data(k,:)                                      =   index_last_peak + subindex_valley_data(k,:) - 1;
%     toc;
    end
    index_valleys_CO2                                           =  index_valleys_CO2(2:length(index_valleys_CO2)-1);

    
% For each dataset, for each peak to the next peak, the data
% points are interpolated into the peak-to-peak coordinates as found in
% CO2. This results in compressed/stretched datapoints between a peak and
% the next peak, matching those peaks of CO2.
    coupled_indices_intensity_peakvalley_data{1,col}           =   [];
    sorted_coupled_indices_intensity_peakvalley_data{1,col}    =   [];
    % Create a sorted list of indices of peaks and valleys as they will
    % appear in the adjusted_intensity_data.
    coupled_indices_intensities_peakvalley                          =   [index_peaks_CO2,intensity_peak_data_absorbance;index_valleys_CO2,intensity_valley_data];
    sorted_coupled_indices_intensity_peakvalley_CO2                   =   sortrows(coupled_indices_intensities_peakvalley);

    for k=1:col;
        coupled_indices_intensity_peakvalley_data{1,k}          =   [index_peak_data(:,k),intensity_peak_data_absorbance(:,k);index_valley_data(:,k),intensity_valley_data(:,k)];
        sorted_coupled_indices_intensity_peakvalley_data{1,k}   =   sortrows(coupled_indices_intensity_peakvalley_data{1,k});
        for m=1:2:length(sorted_coupled_indices_intensity_peakvalley_data{1,k})-2                                                   %Change these twos to ones to also match data valleys to CO2 valleys
            begin_data                                          =   sorted_coupled_indices_intensity_peakvalley_data{1,k}(m,1);
            end_data                                            =   sorted_coupled_indices_intensity_peakvalley_data{1,k}(m+2,1);   %And this two->one
            % Some data points are marked as both a peak and a valley (not
            % sure why), causing a NaN in wavenumber_interval. Following
            % "if" statement circumvents this. 
            if begin_data == end_data
                if m+2 >= length(sorted_coupled_indices_intensity_peakvalley_data{1,k})
                    end_data                                    =   end_data+1;
                else
                    end_data                                        =   sorted_coupled_indices_intensity_peakvalley_data{1,k}(m+2,1);   %Not this two
                end 
            end
            begin_adjust                                                =   sorted_coupled_indices_intensity_peakvalley_CO2(m,1);
            end_adjust                                                  =   sorted_coupled_indices_intensity_peakvalley_CO2(m+2,1); % And this two->one
            wavenumber_interval                                         =   (wavenumber(end_adjust)-wavenumber(begin_adjust))/(end_data-begin_data);
            adjustedCO2_absorbance(begin_adjust:end_adjust,k)  =   interp1(wavenumber(begin_adjust):wavenumber_interval:wavenumber(end_adjust),absorbance(begin_data:end_data,k),wavenumber(begin_adjust:end_adjust),'linear');
        end
    end
    
% Data between the two big CO2 regions gets messed up, so override it here.
% First pick a point between the two regions.
    border_index        =       3400;
    [intensity,index_k_1002]    =   min(abs( ...
        sorted_coupled_indices_intensity_peakvalley_CO2(:,1)-border_index));
    adjustedCO2_absorbance(sorted_coupled_indices_intensity_peakvalley_CO2(index_k_1002-1,1)+1:...
        sorted_coupled_indices_intensity_peakvalley_CO2(index_k_1002+1,1)-1,:)=...
        absorbance(sorted_coupled_indices_intensity_peakvalley_CO2(index_k_1002-1,1)+1:...
        sorted_coupled_indices_intensity_peakvalley_CO2(index_k_1002+1,1)-1,:);


[adjustedH2OCO2_absorbance,index_peaks_H2O,index_H2Opeak_data]= ...
    H2O_calibration(wavenumber,adjustedCO2_absorbance,H2O,index_peaks_CO2(end));
toc

[filename]=WriteCalibratedDataToText(   adjustedH2OCO2_absorbance,...
                                        sampleControlPath,saveTXTsTo,G0S1,...
                                        G0S2,G0S3,G0S4,G1S1,G1S2,G2S1,G2S2);


% close all
% figure;
% absplot     =   2;
% H2O_factor  =   50000;
% hold on
% plot(wavenumber,absorbance(:,absplot), ...
%     wavenumber, adjustedH2OCO2_absorbance(:,absplot),...
% 	wavenumber, CO2*CO2_factor,'',...
% 	wavenumber, H2O*H2O_factor,'r',...
% 	wavenumber(index_peaks_CO2),CO2(index_peaks_CO2)*CO2_factor,'xk',...
% 	wavenumber(index_valleys_CO2),CO2(index_valleys_CO2)*CO2_factor,'xm',...
% 	wavenumber(index_peak_data(:,absplot)),absorbance(index_peak_data(:,absplot),absplot),'+k',...
% 	wavenumber(index_valley_data(:,absplot)),absorbance(index_valley_data(:,absplot),absplot),'+m',...
% 	wavenumber(index_peaks_H2O),H2O(index_peaks_H2O)*H2O_factor,'dk',...
% 	wavenumber(index_H2Opeak_data(:,absplot)),absorbance(index_H2Opeak_data(:,absplot),absplot),'sk')
% hold off
% xlabel('Wavenumber (cm^-^1)')
% ylabel('Absorbance')
% title('Absorbance calibrated to CO_2 and H_2O');
% legend('Absorbance',...
%     'Absorbance adjusted for CO_2 and H_2O',...
%     'CO_2','H_2O',...
%     'CO2 peak identified by the script','CO2 valley identified by the script',...
%     'data peaks matched to CO2','data valleys matched to CO2',...
%     'H2O peak identified by the script','data peaks matched to H2O')
toc
% end