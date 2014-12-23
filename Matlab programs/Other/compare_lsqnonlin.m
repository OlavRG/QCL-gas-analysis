% %   Input
% data_path   =   ['D:\Workspace\Breath Analysis\Measurements\data of 26-8-2014\data no zeros\'];
% list        =   ['G0S1';'G0S2';'G0S3';'G0S4';'G1S1';'G1S2';'G2S1';'G2S2'];
% G0S1        =   importdata([data_path list(1,:) '.txt']); G0S1=G0S1.data;
% G0S2        =   importdata([data_path list(2,:) '.TXT']); G0S2=G0S2.data;
% G0S3        =   importdata([data_path list(3,:) '.TXT']); G0S3=G0S3.data;
% G0S4        =   importdata([data_path list(4,:) '.TXT']); G0S4=G0S4.data;
% G1S1        =   importdata([data_path list(5,:) '.TXT']); G1S1=G1S1.data;
% G1S2        =   importdata([data_path list(6,:) '.TXT']); G1S2=G1S2.data;
% G2S1        =   importdata([data_path list(7,:) '.TXT']); G2S1=G2S1.data;
% G2S2        =   importdata([data_path list(8,:) '.TXT']); G2S2=G2S2.data;

addpath(['D:\Workspace\Breath Analysis\Matlab programs\'...
    'HealthyAsthmaCF analysis\strjoin']);

file_regexp = '.*25T.TXT';
    healthy     =   [G0S1,G0S2];
    asthma      =   [G1S1,G1S2];
    healthy_asthma  =   [healthy, asthma];
    healthy_asthma(isnan(healthy_asthma)) = 0 ;
    options = optimset('Display','off','TolFun',1e-15); % 
    wavenumber= load('D:\Workspace\Breath Analysis\Measurements\Wavenumber.txt');
    wavenumber  =   wavenumber(1:length(G0S1)); 

    INTERACTION_LENGTH = 54.36;
    STANDARD_COMPOUND_PATH  =   'D:\Workspace\Breath Analysis\Compounds\';
    standardCompoundList      =   cellstr(ls(STANDARD_COMPOUND_PATH));
    standardCompoundList2      =   standardCompoundList(3:end,1);
    for k=1:length(standardCompoundList2)
        standardCompoundList{k} = [STANDARD_COMPOUND_PATH standardCompoundList2{k}];
        dirorno(k) = isdir(standardCompoundList{k});
    end
    standardCompoundList = standardCompoundList(dirorno);
    
    standardCompoundDir  =   cell(length(standardCompoundList),1);
    database_compounds_all   =   zeros(length(wavenumber),length(standardCompoundList));
    nStandardCompound=0;
    nEmpty=0;
    % This loop loads and processes the standard compounds to fit the
    % wavenumber length. It also takes into account the possibility of
    % empty folders in STANDARD_COMPOUND_PATH.
    for k=1:length(standardCompoundList)
        standardCompoundDir{k,1}    =   cellstr(ls(... 
            standardCompoundList{k}(1:end)));
        standardFileName   =   strjoin(regexp(standardCompoundDir{k,1},file_regexp,'match'),'');
        if isempty(standardFileName) == 1
            nEmpty=nEmpty+1;
            indexEmptyFolder(nEmpty)     =   k;
        else
            nStandardCompound=nStandardCompound+1;
            [ database_compounds_all(:,nStandardCompound) ] =   load_compound(wavenumber,...
                [standardCompoundList{k}(1:end) '\' standardFileName]);
        end
    end

func_all = @(C) difference_for_least_squares_all( ...
    healthy_asthma(:,1), database_compounds_all*INTERACTION_LENGTH, C);

C_initial_all (1:nStandardCompound,1) = 1;
C_initial_all (5,1) = 1000;     % CO2
C_initial_all (14,1) = 1000;    % Water

lb = zeros(size(C_initial_all)); %lower bound of zero
ub = Inf*ones(size(C_initial_all));

concentrationArray = lsqnonlin(func_all,C_initial_all,lb,ub, options);
