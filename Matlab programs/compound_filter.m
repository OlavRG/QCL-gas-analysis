% This program requires a cell list of compound names as they appear in
% PNNL as input, then filters out the compounds that do not reach
% minimum_absorbance at single wavenumber. The absorbance of each compound
% is computed using an assumed_concentration.

% Can use 'L:\IST\OP\scratch\Olav Grouwstra\Compounds\'; as testing
% directory

clear all
cd('H:\My Documents\GitHub\QCL-gas-analysis');
addpath(genpath('H:\My Documents\GitHub\QCL-gas-analysis\Matlab programs'));
% Find from STANDARD_COMPOUND_PATH all the compounds that at some point have
% a value of minimum_absorbance
    file_regexp     =   '.*_25T?.TXT';
    wavenumber= load('L:\IST\OP\scratch\Olav Grouwstra\Measurements\Wavenumber.txt');
    STANDARD_COMPOUND_PATH  =  'L:\IST\OP\scratch\Adonis\Databases\PNNL Database\Compounds\';
    
% Now let's say there is a concentration of assumed_concentration, then the
% absorbance should be at least minimum_absorbance, which all compounds
% that are present (Absorbance>E-45) in "Absorbance calibrated with
% 6molecules.fig" should get. Higher absorbance and lower concentration
% makes this stricter.
    minimum_absorbance  =   50*10^-7%10^-7; % in base-10 absorbance: I=I0*10^(-A)
    assumed_concentration = 0.001; % in ppmv
    interaction_length = 54.36%1; % in meter

% In case one wishes to filter in a certain wavelength range
%     wavelength_start    = 1020;
%     wavelength_end      = 1100;
     if exist('wavelength_start','var')
         [wavelength_start_index,~]  =   find(wavenumber==wavelength_start);
     else
         wavelength_start_index = 1;
     end
     
     if exist('wavelength_end','var')
         [wavelength_end_index,~]    =   find(wavenumber==wavelength_end);
     else
         wavelength_end_index = length(wavenumber);
     end
     
%     compoundPathFileList      =   cellstr(ls(STANDARD_COMPOUND_PATH));
%     compoundPathFileList      =   compoundPathFileList(3:end,1);
    load('H:\My Documents\GitHub\QCL-gas-analysis\Breath molecules_2.mat')
    compoundPathFileList      =   match_Volatinome_PNNL(:,2);

    
    standardCompoundDir  =   cell(length(compoundPathFileList),1);
    absorptivity_all_compounds   =   zeros(length(wavenumber),length(compoundPathFileList));
    nStandardCompound=0;
    nEmpty1=0;
    indexEmptyFolder1 = zeros(length(compoundPathFileList),1);
    % This loop loads and processes the standard compounds to fit the
    % wavenumber length. It also takes into account the possibility of
    % folders with no file matching the regex pattern in
    % STANDARD_COMPOUND_PATH.
    for k=1:length(compoundPathFileList)
        standardCompoundDir{k,1}    =   cellstr(ls(fullfile(STANDARD_COMPOUND_PATH, ... 
            compoundPathFileList{k}(1:end))));
        standardFileName   =   strjoin(regexpi(standardCompoundDir{k,1},file_regexp,'match'),'');
        if isempty(standardFileName) == 1
            nEmpty1=nEmpty1+1;
            indexEmptyFolder1(nEmpty1)     =   k;
        else
            nStandardCompound=nStandardCompound+1;
            [ absorptivity_all_compounds(:,nStandardCompound) ] =   load_compound(wavenumber,...
                [STANDARD_COMPOUND_PATH compoundPathFileList{k}(1:end) '\' standardFileName]);
        end
    end

    % Resize absorptivity_all_compounds to account for folders without a
    % file matching the regex. 
    absorptivity_all_compounds  =   absorptivity_all_compounds(:,1:nStandardCompound);
    standardCompoundList = compoundPathFileList;
    indexEmptyFolder1(indexEmptyFolder1==0) = [];
    if exist('indexEmptyFolder1','var') 
        standardCompoundList(indexEmptyFolder1)=[];
    end
    [row_standard_compounds, nStandardCompound] = size(absorptivity_all_compounds);
    
    database_compounds_all_absorbance = absorptivity_all_compounds*interaction_length*assumed_concentration;
    absorptivity_all_compounds_filtered = zeros(row_standard_compounds,nStandardCompound);
    nPassedCompound = 0;
    nCutCompound = 0;
    indexEmptyFolder2 = zeros(length(standardCompoundList),1);
    for ii=1:nStandardCompound
        if max(database_compounds_all_absorbance(wavelength_start_index:wavelength_end_index,ii)) >= minimum_absorbance
            nPassedCompound = nPassedCompound + 1;
            absorptivity_all_compounds_filtered(:,nPassedCompound) = absorptivity_all_compounds(:,ii);
        else
            nCutCompound = nCutCompound + 1;
            indexEmptyFolder2(nCutCompound)     =   ii;
        end
    end
    
    absorptivity_all_compounds_filtered = absorptivity_all_compounds_filtered(:,1:nPassedCompound);
    filteredCompoundList = standardCompoundList;
    indexEmptyFolder2(indexEmptyFolder2==0) = [];
    if exist('indexEmptyFolder2','var') 
        filteredCompoundList(indexEmptyFolder2)=[];
    end
% save compounds_twice_filtered.mat absorptivity_all_compounds_filtered filteredCompoundList wavenumber

% Also write the variables to .txt files for processing outside of MATLAB
% dlmwrite('absorptivity_all_compounds_filtered.txt',absorptivity_all_compounds_filtered,'precision','%.15E','delimiter','\t');
% fileID1=fopen('filteredCompoundList.txt','w')
% [rowCompList, colCompList] = size(filteredCompoundList);
% for row=1:rowCompList
%     fprintf(fileID1,'%s\n', filteredCompoundList{row,:});
% end
% fclose(fileID1);