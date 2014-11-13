

%Dependencies: 
%   difference_for_least_squares_all.m       
%   compound_region.mat
%       compound_path
%       file_regexp
%       file_extension
%       p_region
%       compound_region
%   Absorbances.mat   absorbance measurements of all samples (output by CO2
%   calibration script)
%   Standard compound files and path
%   compound files of compounds in compound_region

% %   Input
% data_path   =   ['L:\IST\OP\scratch\Olav Grouwstra\Measurements\'...
%                 'CO2H2O calibrated data\data of 26-8-2014\data no zeros\'];
% list        =   ['G0S1';'G0S2';'G0S3';'G0S4';'G1S1';'G1S2';'G2S1';'G2S2'];
% G0S1        =   importdata([data_path list(1,:) '.TXT']); G0S1=G0S1.data;
% G0S2        =   importdata([data_path list(2,:) '.TXT']); G0S2=G0S2.data;
% G0S3        =   importdata([data_path list(3,:) '.TXT']); G0S3=G0S3.data;
% G0S4        =   importdata([data_path list(4,:) '.TXT']); G0S4=G0S4.data;
% G1S1        =   importdata([data_path list(5,:) '.TXT']); G1S1=G1S1.data;
% G1S2        =   importdata([data_path list(6,:) '.TXT']); G1S2=G1S2.data;
% G2S1        =   importdata([data_path list(7,:) '.TXT']); G2S1=G2S1.data;
% G2S2        =   importdata([data_path list(8,:) '.TXT']); G2S2=G2S2.data;


%% v2.1 Ammonia (NH3), methane (CH4) and ethanol (ETOH) added (Adonis)
%   v2.2 Loading/processing of standard compounds is done through
%   load_compound.m (Olav)

%%

timeId1=tic;
clear('database_compounds_all','func_all','lb','ub','C_initial_all','concentrations');


addpath(['L:\IST\OP\scratch\Olav Grouwstra\Matlab programs\'...
    'HealthyAsthmaCF analysis\strjoin']);
    %Added so downloaded strjoin is loaded over strjoin of matlab's toolbox. 

% load Healthy_Asthma_Raw_Calibrated.mat %Contains multiple variables
% load Absorbances_raw_20140826_Calibrated_NoZeros.mat;
    %Also contains Wavenumber vector
    healthy     =   [G0S1,G0S2];
    asthma      =   [G1S1,G1S2];
    healthy_asthma  =   [healthy, asthma];
    healthy_asthma(isnan(healthy_asthma)) = 0 ;
    options = optimset('Display','off','TolFun',1e-15); % 
    wavenumber= load('L:\IST\OP\scratch\Olav Grouwstra\Measurements\Wavenumber.txt');
    wavenumber  =   wavenumber(1:length(G0S1)); 
    
% Establish some sizes for later use
    [row_healthy,samples_healthy] =   size(healthy);
    [row_asthma,samples_asthma] =   size(asthma);
    [row_p_region,col_p_region] =   size(p_region);
    [row_compound_region,col_compound_region]   =   size(compound_region);

% This counts the compounds that are significantly present in
% region_pres or more regions, so preallocation can be done for next loop.
    region_pres=1; % Max value is col_p_region;
    nExtraCompound=sum(sum(cell2mat(...
        compound_region(2:row_compound_region,2:1+col_p_region)),2) >= region_pres);

% Prepare H2O and CO2 acetone and ammonia
    INTERACTION_LENGTH = 54.36;
    STANDARD_COMPOUND_PATH  =   'L:\IST\OP\scratch\Olav Grouwstra\Compounds\';
    standardCompoundList      =   cellstr(ls(STANDARD_COMPOUND_PATH));
    standardCompoundList      =   standardCompoundList(3:end,1);

    standardCompoundDir  =   cell(length(standardCompoundList),1);
    database_compounds_all   =   zeros(length(wavenumber),length(standardCompoundList));
    nStandardCompound=0;
    nEmpty=0;
    % This loop loads and processes the standard compounds to fit the
    % wavenumber length. It also takes into account the possibility of
    % empty folders in STANDARD_COMPOUND_PATH.
    for k=1:length(standardCompoundList)
        standardCompoundDir{k,1}    =   cellstr(ls(fullfile(STANDARD_COMPOUND_PATH, ... 
            standardCompoundList{k}(1:end))));
        standardFileName   =   strjoin(regexp(standardCompoundDir{k,1},file_regexp,'match'),'');
        if isempty(standardFileName) == 1
            nEmpty=nEmpty+1;
            indexEmptyFolder(nEmpty)     =   k;
        else
            nStandardCompound=nStandardCompound+1;
            [ database_compounds_all(:,nStandardCompound) ] =   load_compound(wavenumber,...
                [STANDARD_COMPOUND_PATH standardCompoundList{k}(1:end) '\' standardFileName]);
        end
    end

    % Resize database_compounds_all to account for empty folders in 
    % length(standardCompoundList) to which it was calibrated
    database_compounds_all  =   database_compounds_all(:,1:nStandardCompound);
    if exist('indexEmptyFolder','var') 
        standardCompoundList(indexEmptyFolder)=[];
    end
    [row_standard_compounds, nStandardCompound] = size(database_compounds_all);


    % Preallocate for coming loop
    compound_files    =   cell(nExtraCompound,1);
    database_compounds_all(:,nStandardCompound+1:nStandardCompound+nExtraCompound)...
        =   zeros(length(wavenumber),nExtraCompound);
    
    toc(timeId1);
j=0;
for k=2:4%row_compound_region
    if sum([compound_region{k,2:1+col_p_region}])     >=   region_pres
        j=j+1;
        compound_files{j,1}   =   cellstr(ls(fullfile(compound_path, ... 
            compound_region{k,1}(1:end))));
        standardFileName   =   strjoin(regexp(compound_files{j,1},file_regexp,'match'),'');
        database_compounds_all(:,nStandardCompound+j) = load_compound ...
    ( wavenumber, [compound_path compound_region{k,1}(1:end) '\' standardFileName]);
        
    time1=toc(timeId1);
    disp([num2str(time1) 's ' num2str(k) ' out of ' num2str(row_compound_region)]) 
    else
    end
    
end

%allocate memory for the concentrations
    concentrationArray = zeros(samples_healthy+samples_asthma,nExtraCompound+nStandardCompound);

    %Initial concentration guess in ppm.
    C_initial_all (1:nStandardCompound,1) = 1;
    C_initial_all (5,1) = 1000;     % CO2
    C_initial_all (14,1) = 1000;    % Water
    C_initial_all (nStandardCompound+1:nExtraCompound+nStandardCompound,1) = 0.01;
 
    lb = zeros(size(C_initial_all)); %lower bound of zero
    ub = Inf*ones(size(C_initial_all));

for i=1:samples_healthy+samples_asthma;
    
    func_all = @(C) difference_for_least_squares_all( ...
        healthy_asthma(:,i), database_compounds_all*INTERACTION_LENGTH, C);

    concentrationArray(i,:) = lsqnonlin(func_all,C_initial_all,lb,ub, options);
    time1=toc(timeId1);
    disp([num2str(time1) 's ' num2str(i) ' out of ' num2str(samples_healthy+samples_asthma)]) 
end

compoundListTotal = [standardCompoundList;compound_region(2:end,1)];
% dlmwrite('concentrationArray.txt',concentrationArray,'precision','%1.15g')


% figure;
% dg = [0 0.5 0];
% sample=1;
% absoverlaid=plot(wavenumber,healthy_asthma(:,sample),...
%     wavenumber,database_compounds_all(:,1)*interaction_length*mean_concentration_raw_healthy_asthma{sample+1,1+1},...%,'Color',dg)%,...
%     wavenumber,database_compounds_all(:,2)*interaction_length*mean_concentration_raw_healthy_asthma{sample+1,2+1},...
%     wavenumber,database_compounds_all(:,3)*interaction_length*mean_concentration_raw_healthy_asthma{sample+1,3+1},...
%     wavenumber,database_compounds_all(:,4)*interaction_length*mean_concentration_raw_healthy_asthma{sample+1,4+1},...
%     wavenumber,database_compounds_all(:,5)*interaction_length*mean_concentration_raw_healthy_asthma{sample+1,5+1},...
%     wavenumber,database_compounds_all(:,6)*interaction_length*mean_concentration_raw_healthy_asthma{sample+1,6+1})%,'Color',[0.5 0.5 0.5])
% xlabel('Wavenumber (cm^-^1)')
% ylabel('Absorbance')
% title('Absorbance spectrum with molecule spectra overlaid');
% legend(absoverlaid, 'Absorbance',...
%     'CO_2',...
%     'H_2O',...
%     'Acetone',...
%     'Ammonia',...
%     'Methane',...
%     'Ethanol')
% 
% %Separate healthy and asthmatic data and get 
% concentrations_healthy = concentrationArray(1:samples_healthy,:);
% concentrations_asthma = concentrationArray(samples_healthy+1:samples_healthy+samples_asthma,:);
% 
% mean_concentrations_healthy = nanmean(concentrations_healthy);
% mean_concentrations_asthma = nanmean(concentrations_asthma);
% std_concentrations_healthy = std(concentrations_healthy);
% std_concentrations_asthma = std(concentrations_asthma);
%         
% mean_concentration_filtered{1,1}    =   ['Mean concentrations (ppm)'];
% mean_concentration_filtered{2,1}    =   ['Healthy mean'];
% mean_concentration_filtered{3,1}    =   ['Healthy st. dev.'];
% mean_concentration_filtered{4,1}    =   ['Asthma mean'];
% mean_concentration_filtered{5,1}    =   ['Asthma st. dev.'];
% 
% mean_concentration_filtered{1,2}    =   ['CO2'];
% mean_concentration_filtered{1,3}    =   ['H2O'];
% mean_concentration_filtered{1,4}    =   ['ACETONE'];
% mean_concentration_filtered{1,5}    =   ['NH3'];
% mean_concentration_filtered{1,6}    =   ['CH4'];
% mean_concentration_filtered{1,7}    =   ['ETHANOL'];
% 
% comp_number = num_compounds+nStandardCompound;
% mean_concentration_filtered(1,nStandardCompound+2:comp_number+1)  =   sample_concentration(1,:);
% %mean_concentration_filtered(1,2:comp_number+1)  =   num2cell(sample_concentration(1,:));
% mean_concentration_filtered(2,2:comp_number+1)    =   num2cell(mean_concentrations_healthy);
% mean_concentration_filtered(3,2:comp_number+1)    =   num2cell(std_concentrations_healthy);
% mean_concentration_filtered(4,2:comp_number+1)    =   num2cell(mean_concentrations_asthma);
% mean_concentration_filtered(5,2:comp_number+1)    =   num2cell(std_concentrations_asthma);
% 
% mean_concentration_raw_healthy_asthma = mean_concentration_filtered;
% 
% save Healthy_Asthma_Raw_Concentrations_Calibrated.mat mean_concentration_raw_healthy_asthma concentrations_asthma...
%     concentrations_healthy





