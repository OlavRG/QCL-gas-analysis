function [ compound5 ] = load_compound ...
    ( wavenumber, string_compound)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     compound1=importdata(string_compound);
    fid = fopen(string_compound);
    compound1 = textscan(fid, '%f %f');
    fclose(fid);

% Convert HITRAN's base-e absorption coeff alpha to base-10 absorbance A
%     if isstruct(compound1)   %used for identifying HITRAN TXT's when using MATLAB's load function
%         compound2   =   -log10(exp(-compound1.data(:,2)*10^-4));
%         compound3     =   [compound1.data(:,1),compound2];
%     else %PNNL TXT's have class double
        compound3=cell2mat(compound1) ;
%     end

    [~,compound_start_index]    =   ...
        min(abs(compound3(:,1)-wavenumber(1)));
    [~,compound_end_index]    =   ...
        min(abs(compound3(:,1)-wavenumber(end)));

% This if is used to flip the compound data vectors, which
% tend to have the highest k in (1,1), and the lowest in (end) for PNNL
% data.
    if compound_start_index > compound_end_index ;
        compound4   =   flipdim(compound3,1);
        lower_compound_index     =   length(compound4)-compound_start_index+1-1;
        upper_compound_index     =   length(compound4)-compound_end_index+1+1;
    else
        compound4   =   compound3;
        lower_compound_index     =   compound_start_index-1;
        upper_compound_index     =   compound_end_index+1;
    end

    compound5     =   interp1( ...
    compound4(lower_compound_index:upper_compound_index,1), ...
    compound4(lower_compound_index:upper_compound_index,2), ...
    wavenumber);


end
