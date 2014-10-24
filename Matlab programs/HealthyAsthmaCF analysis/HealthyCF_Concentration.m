%Dependencies: 
%   difference_for_least_squares_triple.m       
%   compound_region.mat
%       compound_path
%       file_regexp
%       file_extension
%       p_region
%       compound_region
%   raw_matlab.mat                              (k-absorbance measurements)
%   CO2_25T.TXT
%   H2O_25T.TXT
%   compound files of compounds in compound_region
tic

load CF_compound_region.mat %Contains multiple variables
load 'H:\MEP\Measurements\Final Data TUDelft\raw_matlab.mat';
    %Also contains Wavenumber vector
    healthy     =   [G0S1,G0S2];
    CF      =   [G2S1,G2S2];
    healthy_CF  =   [healthy, CF];
    options = optimset('Display','off','TolFun',1e-15); % 

% Establish some sizes for later use
    [row_healthy,samples_healthy] =   size(healthy);
    [row_CF,samples_CF] =   size(CF);
    [row_p_region,col_p_region] =   size(p_region);
    [row_compound_region,col_compound_region]   =   size(compound_region);

% This loop counts the compounds that are significantly present in all
% regions, so preallocation can be done for next loop.
    count1=0;
    for k=2:row_compound_region
        if sum([compound_region{k,2:1+col_p_region}])     == ...
                col_p_region-1
            count1=count1+1;
        else
        end
    end


% Prepare H2O and CO2 and acetone
    interaction_length = 54.36;
    load(['L:\IST\OP\scratch\Adonis\Databases\PNNL Database\Compounds\' ...
        'Carbon_dioxide\CO2_25T.TXT'])
    load(['L:\IST\OP\scratch\Olav Grouwstra\H2O_HITRAN_base10_1ppm_1meter'])
    load(['L:\IST\OP\scratch\Adonis\Databases\PNNL Database\Compounds\'...
        'Acetone\ACETONE_25T.TXT'])
    H2O_25T     =   [Wavenumber,H2O_base10_1ppm_1meter];
    [wave_start,CO2_25T_start_index]    =   ...
        min(abs(CO2_25T(:,1)-Wavenumber(1)));
    [wave_end,CO2_25T_end_index]    =   ...
        min(abs(CO2_25T(:,1)-Wavenumber(end)));

    [wave_start,H2O_25T_start_index]    =   ...
        min(abs(H2O_25T(:,1)-Wavenumber(1)));
    [wave_end,H2O_25T_end_index]    =   ...
        min(abs(H2O_25T(:,1)-Wavenumber(end)));

    [wave_start,ACETONE_25T_start_index]    =   ...
        min(abs(ACETONE_25T(:,1)-Wavenumber(1)));
    [wave_end,ACETONE_25T_end_index]    =   ...
        min(abs(ACETONE_25T(:,1)-Wavenumber(end)));
    
% These and further if's are used to flip the compound data vectors,
% which tend to have the highest k in (1,1), end the lowest in (end).
    if CO2_25T_start_index > CO2_25T_end_index ;
        CO2_25T   =   flipdim(CO2_25T,1);
        lower_CO2_25T_index     =   length(CO2_25T)-CO2_25T_start_index+1-1;
        upper_CO2_25T_index     =   length(CO2_25T)-CO2_25T_end_index+1+1;
    else
        lower_CO2_25T_index     =   CO2_25T_start_index-1;
        upper_CO2_25T_index     =   CO2_25T_end_index+1;
    end
    
    if H2O_25T_start_index > H2O_25T_end_index ;
        H2O_25T   =   flipdim(H2O_25T,1);
        lower_H2O_25T_index     =   length(H2O_25T)-H2O_25T_start_index+1-1;
        upper_H2O_25T_index     =   length(H2O_25T)-H2O_25T_end_index+1+1;
    else
        lower_H2O_25T_index     =   H2O_25T_start_index;%-1;
        upper_H2O_25T_index     =   H2O_25T_end_index;%+1;
    end

    if ACETONE_25T_start_index > ACETONE_25T_end_index ;
        ACETONE_25T   =   flipdim(ACETONE_25T,1);
        lower_ACETONE_25T_index     =   length(ACETONE_25T)-ACETONE_25T_start_index+1-1;
        upper_ACETONE_25T_index     =   length(ACETONE_25T)-ACETONE_25T_end_index+1+1;
    else
        lower_ACETONE_25T_index     =   ACETONE_25T_start_index;%-1;
        upper_ACETONE_25T_index     =   ACETONE_25T_end_index;%+1;
    end

    CO2_25T_r     =   interp1( ...
        CO2_25T(lower_CO2_25T_index:upper_CO2_25T_index,1), ...
        CO2_25T(lower_CO2_25T_index:upper_CO2_25T_index,2), ...
        Wavenumber);
    
    H2O_25T_r     =   interp1( ...
        H2O_25T(lower_H2O_25T_index:upper_H2O_25T_index,1), ...
        H2O_25T(lower_H2O_25T_index:upper_H2O_25T_index,2), ...
        Wavenumber);

    ACETONE_25T_r     =   interp1( ...
        ACETONE_25T(lower_ACETONE_25T_index:upper_ACETONE_25T_index,1), ...
        ACETONE_25T(lower_ACETONE_25T_index:upper_ACETONE_25T_index,2), ...
        Wavenumber);

% Preallocate for coming loop
    compound_files    =   cell(count1,1);
    sample_concentration    =   cell(1+samples_healthy+samples_CF,count1);
    CO2_concentration       =   cell(1+samples_healthy+samples_CF,count1);
    H2O_concentration       =   cell(1+samples_healthy+samples_CF,count1);
    acetone_concentration       =   cell(1+samples_healthy+samples_CF,count1);

j=0;
for k=2:row_compound_region
    if sum([compound_region{k,2:1+col_p_region}])     ==    col_p_region-1 
        j=j+1;
        sample_concentration{1,j}   =   compound_region{k,1}(1:end);
        H2O_concentration{1,j}   =   compound_region{k,1}(1:end);
        CO2_concentration{1,j}   =   compound_region{k,1}(1:end);
        acetone_concentration{1,j}   =   compound_region{k,1}(1:end);
        compound_files{j,1}   =   cellstr(ls(fullfile(compound_path, ... 
            compound_region{k,1}(1:end))));
        file_name   =   strjoin(regexp(compound_files{j,1},file_regexp,'match'),'');
        load([compound_path compound_region{k,1}(1:end) '\' file_name]);
            % For speed: comment "load" iff files already loaded
        compound_copy    =   eval(regexprep(file_name,file_extension,''));
        % Find overlapping region of database compound (compound_copy) and 
        % measured compounds (Wavenumber), and resample the database
        % compound.
        [wave_start,compound_start_index]    =   ...
            min(abs(compound_copy(:,1)-Wavenumber(1)));
        [wave_end,compound_end_index]    =   ...
            min(abs(compound_copy(:,1)-Wavenumber(end)));
        if compound_start_index > compound_end_index;
            compound_copy   =   flipdim(compound_copy,1);
            lower_k_index     =   length(compound_copy)-compound_start_index+1-1;
            upper_k_index     =   length(compound_copy)-compound_end_index+1+1;
        else
            lower_k_index     =   compound_start_index-1;
            upper_k_index     =   compound_end_index+1;
        end
%         compound_copy_r     =   resample( ...
%             compound_copy(lower_k_index:upper_k_index,2), ...
%             length(Wavenumber), ...
%             upper_k_index-lower_k_index+1);
        compound_copy_r     =   interp1( ...
            compound_copy(lower_k_index:upper_k_index,1), ...
            compound_copy(lower_k_index:upper_k_index,2), ...
            Wavenumber);
        for i=1:samples_healthy+samples_CF;
            
            func = @(C) difference_for_least_squares_triple( ...
                healthy_CF(:,i), compound_copy_r*interaction_length, ...
                CO2_25T_r*interaction_length, ...
                H2O_25T_r*interaction_length, ...
                ACETONE_25T_r*interaction_length, C );
            C_initial = [0.001 1000 1000 1]; %Initial concentration guess in ppm.
            temp = lsqnonlin(func,C_initial,[],[], ...
                options);
            sample_concentration{1+i,j}     =   temp(1);
            CO2_concentration{1+i,j}        =   temp(2);
            H2O_concentration{1+i,j}        =   temp(3);
            acetone_concentration{1+i,j}        =   temp(4);
        end
toc
    else
    end
end

mean_concentration{1,1}    =   ['Mean concentrations (ppm)'];
mean_concentration{2,1}    =   ['Healthy mean'];
mean_concentration{3,1}    =   ['Healthy st. dev.'];
mean_concentration{4,1}    =   ['CF mean'];
mean_concentration{5,1}    =   ['CF st. dev.'];
    
mean_concentration(1,2:count1+1)  =   sample_concentration(1,:);

mean_concentration(2,2:count1+1)    =   num2cell(mean(cell2mat(sample_concentration ...
    (2:samples_healthy+1,:))));
mean_concentration(3,2:count1+1)    =   num2cell(std(cell2mat(sample_concentration ...
    (2:samples_healthy+1,:))));
mean_concentration(4,2:count1+1)    =   num2cell(mean(cell2mat(sample_concentration ...
    (samples_healthy+2:samples_CF+samples_healthy+1,:))));
mean_concentration(5,2:count1+1)    =   num2cell(std(cell2mat(sample_concentration ...
    (samples_healthy+2:samples_CF+samples_healthy+1,:))));

plot(Wavenumber,healthy_CF(:,i), ...
    Wavenumber,compound_copy_r*temp(1)*interaction_length, ...
    Wavenumber,CO2_25T_r*temp(2)*interaction_length, ...
    Wavenumber,H2O_25T_r*temp(3)*interaction_length, ...
    Wavenumber,ACETONE_25T_r*temp(4)*interaction_length, ...
    Wavenumber,smooth(healthy_CF(:,i),100))
%     Wavenumber,(compound_copy_r*temp(1)+CO2_25T_r*temp(2)+ ...
%     H2O_25T_r*temp(3))*interaction_length)
legend('Measurement','m-Cresol','CO2','H2O','Acetone','smooth')%,'Sum of compounds')



