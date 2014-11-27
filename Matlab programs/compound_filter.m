% This program requires the folder of saved PNNL compounds as input, then
% filters out the compounds that at do not reach minimum_absorbance at
% single wavenumber. The absorbance of each compound is computed using an
% assumed_concentration.

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
    minimum_absorbance  =   10^-7; % in base-10 absorbance: I=I0*10^(-A)
    assumed_concentration = 0.001; % in ppmv
    INTERACTION_LENGTH = 1; % in meter

%     compoundPathFileList      =   cellstr(ls(STANDARD_COMPOUND_PATH));
%     compoundPathFileList      =   compoundPathFileList(3:end,1);
    load('H:\My Documents\GitHub\QCL-gas-analysis\Breath molecules_2.mat')
    compoundPathFileList      =   match_Volatinome_PNNL(:,2);

    
    standardCompoundDir  =   cell(length(compoundPathFileList),1);
    database_compounds_all   =   zeros(length(wavenumber),length(compoundPathFileList));
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
            [ database_compounds_all(:,nStandardCompound) ] =   load_compound(wavenumber,...
                [STANDARD_COMPOUND_PATH compoundPathFileList{k}(1:end) '\' standardFileName]);
        end
    end

    % Resize database_compounds_all to account for empty folders in 
    % length(standardCompoundList) to which it was calibrated
    database_compounds_all  =   database_compounds_all(:,1:nStandardCompound);
    standardCompoundList = compoundPathFileList;
    indexEmptyFolder1(indexEmptyFolder1==0) = [];
    if exist('indexEmptyFolder1','var') 
        standardCompoundList(indexEmptyFolder1)=[];
    end
    [row_standard_compounds, nStandardCompound] = size(database_compounds_all);
    
    database_compounds_all2 = database_compounds_all*INTERACTION_LENGTH*assumed_concentration;
    database_compounds_all3 = zeros(row_standard_compounds,nStandardCompound);
    nPassedCompound = 0;
    nCutCompound = 0;
    indexEmptyFolder2 = zeros(length(standardCompoundList),1);
    for ii=1:nStandardCompound
        if max(database_compounds_all2(:,ii)) >= minimum_absorbance
            nPassedCompound = nPassedCompound + 1;
            database_compounds_all3(:,nPassedCompound) = database_compounds_all2(:,ii);
        else
            nCutCompound = nCutCompound + 1;
            indexEmptyFolder2(nCutCompound)     =   ii;
        end
    end
    
    filteredCompoundList = standardCompoundList;
    indexEmptyFolder2(indexEmptyFolder2==0) = [];
    if exist('indexEmptyFolder2','var') 
        filteredCompoundList(indexEmptyFolder2)=[];
        database_compounds_all3(:,indexEmptyFolder2)=[];
    end
save compounds_twice_filtered.mat database_compounds_all3 filteredCompoundList
    