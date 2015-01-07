% Use:
%       Finds regions of significant difference between healthy and
%       asthmatic patients based purely on p-values. Then checks which
%       compounds have a minimal of min_intensity in these regions, and
%       lists them in compound_region. Also saves a list of variables to
%       compound_region.mat for use in HealthyAsthma_concentrations.m.

% Input:
%       All compounds
%       p-values of healthy vs asthma
%       Wavenumber vector

% Output:
%       compound_region.mat


% Dependencies:
% strjoin (CAUTION: not the standard matlab strjoin.)

addpath(genpath('D:\Workspace\Breath Analysis\Matlab programs\HealthyAsthmaCF analysis\strjoin'))

tic;
close all
clear all

% Declare variables
    max_p           =   0.05;
    min_p_region    =   10;
    min_intensity   =   10^(-6);
    interaction_length = 54.36; % in meters
    compound_path   =   'W:\';%['L:\IST\OP\scratch\Adonis\Databases\' ... 
        %'PNNL Database\Compounds\'];
    file_regexp     =   '.*_25T?.TXT';
    file_extension  =   '.TXT';
% Set indices of folders of compounds to use
    compound_list   =   cellstr(ls(compound_path));
    compound_list   =   compound_list(3:end,1);
    lower_compound_index    =   1;
    upper_compound_index    =   length(compound_list);%length(compound_list);

        
% load and analyse dataset, load p values and corresponding wavenumber vector
%     path_to_raw_matlab  =   ['H:\MEP\Measurements\Final Data TUDelft\' ... 
%         'raw_matlab.mat'];
%     [healthy_mean_std, asthma_mean_std, cysticF_mean_std]  =    ...
%         Mean_std_analysis(path_to_raw_matlab);
%     pvalues=    importdata(['L:\IST\OP\scratch\Olav Grouwstra\Measurements\'...
%         'P-values of CO2H2O calibrated data\'...
%         'raw.quantile.healthy.vs.asthma.alldata.26-8-2014.no zeros.csv'],'\t',1);
    pvalues=    importdata(['D:\Workspace\Breath Analysis\'...
        'Measurements\P-values of CO2H2O calibrated data\'...
        'raw.MAD.healthy.vs.CF.alldata.csv'],',',1);

    wavenumber= load('D:\Workspace\Breath Analysis\Measurements\Wavenumber.txt');
    
% Determine indices of wavenumber regions of at least min_p_region 
% consecutive p-values < max_p
    [ind]                   =   find(pvalues.data(:,5) <= max_p);
    length_ind              =   length(ind);
    ind_filtered            =   zeros(length_ind,1);
    for i=2:length_ind-min_p_region+1               % Make work for i=1:len...
        if ind(i+min_p_region-1)-ind(i)== min_p_region-1   
            ind_filtered(i:i+min_p_region-1,1)  =   ... 
                ind(i:i+min_p_region-1,1);
        end
    end
    ind_filtered=ind_filtered(ind_filtered~=0);  % Remove zeros
    mIndFiltered=length(ind_filtered);
    
% Establish the begin and end indices of wavenumber regions in matrices
    ind_filtered_shift(1:mIndFiltered-1,1)    =   ind_filtered(2:mIndFiltered);
    shift_filtered  =     ind_filtered_shift-ind_filtered(1:mIndFiltered-1);
    ind_filtered_begin              =   find( ...
        shift_filtered > 1)+1;
    ind_filtered_begin              =   [1; ind_filtered_begin];
    ind_filtered_end                =   [ind_filtered_begin(2:end)-1; mIndFiltered];

    region_amount   =   length(ind_filtered_begin);
    
% Put wavenumber regions of at least min_p_region consecutive 
% p-values < max_p in cells
    p_region    =   cell(1,length(ind_filtered_begin));
    for k=1:length(ind_filtered_begin)
        p_region{1,k}   =   wavenumber(ind_filtered( ...
            ind_filtered_begin(k):ind_filtered_end(k)));
    end

% Allocate table to display what compounds are in what region, compound_region
    compound_amount         =   (upper_compound_index-lower_compound_index)+1;
    compound_region         =   num2cell(zeros(compound_amount+1, ...
        length(ind_filtered_begin)+1));
    compound_region{1,1}    =   ['Compounds\Wavenumber cm-1'];

% Get list of the folders in which compound files are located
    index_start     =   zeros(compound_amount,length(ind_filtered_begin));
    index_end     =   zeros(compound_amount,length(ind_filtered_begin));
    
% For each compound:
    % Put list of files of the compound in compound_dir, find the right
    % file and load it. To access variable with string name, copy it to
    % temporary compound_copy. Now for each region find the start and end
    % index for the particular compound, and if all intensities within the
    % region > min_intensity, then set it to one in the compound_region
    % table
    
% To be increased in speed by saving compounds to a single binary matfile 
% initially, and then loading once for actual runs:
% save myfile.mat <var_a> <var_b> ...
    compound_dir    =   cell(compound_amount,1);
    for k=lower_compound_index:upper_compound_index
        compound_dir{k,1}    =   cellstr(ls(fullfile(compound_path, ... 
            compound_list{k}(1:end))));
        file_name   =   strjoin(regexp(compound_dir{k,1},file_regexp,'match'),''); 
        if isempty(file_name) == 1      %Above line FAILS when using MATLAB's strjoin
        
        else
%             load([compound_path compound_list{k}(1:end) '\' file_name]);
            fid = fopen([compound_path compound_list{k}(1:end) '\' file_name]);
            compound_copy = textscan(fid, '%f %f');
            fclose(fid);
            compound_copy = cell2mat(compound_copy);
%             file_name_wo_ext = regexprep(file_name,file_extension,'');
%             var_name = regexprep(file_name_wo_ext,'(^[0-9]\w+)','X$1');
%             compound_copy    =   eval(regexprep(var_name,file_extension,''));
            compound_region{k+1,1}  =   compound_list{k};
            for j=1:length(ind_filtered_begin)
                compound_region{1,j+1}      =   [num2str(p_region{j}(1)) '-' ... 
                    num2str(p_region{j}(end))];
                [k_dif_start,index_start(k,j)]   =   min(abs(p_region{j}(1)-compound_copy(:,1)));
                [k_dif_end,index_end(k,j)]       =   min(abs(p_region{j}(end)-compound_copy(:,1)));
                if mean(compound_copy(index_end(k,j):index_start(k,j),2) > min_intensity) == 1
                    compound_region{k+1,j+1}    =   1;
                end
            end
        end
        toc
    end

    
% This while loop removes all compound entries that do not sufficiently
% match any region (all zeros in compound_region).
    k=1;
    [rows_compound_region,cols_compound_region]    =   size(compound_region);
    while k<(rows_compound_region)
        if sum([compound_region{k+1,2:length(ind_filtered_begin)+1}])==0
            compound_region(k+1,:)=[];
        else
            k=k+1;
        end
    [rows_compound_region,cols_compound_region]    =   size(compound_region);
    end
beep
% Save variables for Concentration_determination.m
% save raw.MAD.healthy.vs.CF.alldata.alldata_compound_region_min_int_E-6.mat compound_region p_region compound_path file_regexp ...
%     file_extension
%     
toc