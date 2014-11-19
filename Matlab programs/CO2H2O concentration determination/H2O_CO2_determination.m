

%Dependencies: 
%   difference_for_least_squares_all.m       
%   CO2_HITRAN.TXT
%   H2O_HITRAN.TXT
%   absorbance_mean_adjusted_data
%   absorbance_mean_data_old
%tic

% Made by Olav for running after Raw_data_to_absorbance.m to plot H2O and
% CO2 on the calibrated absorbance.


clear('database_compounds_all','func_all','lb','ub','C_initial_all','concentrations');

addpath(['L:\IST\OP\scratch\Olav Grouwstra\Matlab programs\strjoin']);
    %Added so downloaded strjoin is loaded over strjoin of matlab's toolbox. 

    %Also contains wavenumber_data vector
    healthy     =   [absorbance_mean_adjusted_data];
    asthma      =   [absorbance_mean_data_old];
    healthy_asthma  =   [healthy, asthma];
    healthy_asthma(isnan(healthy_asthma)) = 0 ;
    options = optimset('Display','off','TolFun',1e-15); % 
   

% Establish some sizes for later use
    [row_healthy,samples_healthy] =   size(healthy);
    [row_asthma,samples_asthma] =   size(asthma);

% Prepare H2O and CO2 acetone and ammonia
    interaction_length = 54.36;
    CO2_HITRAN = importdata(['L:\IST\OP\scratch\Olav Grouwstra\Compounds\' ...
        'CO2_HITRAN.TXT']);   
    H2O_HITRAN = importdata(['L:\IST\OP\scratch\Olav Grouwstra\Compounds\' ...
        'H2O_HITRAN.TXT']);   

% Convert HITRAN's base-e absorption coeff alpha to base-10 absorbance A
    CO2_HITRAN.data(:,2)=-log10(exp(-CO2_HITRAN.data(:,2)*10^-4));
    CO2_25T     =   [CO2_HITRAN.data(:,1),CO2_HITRAN.data(:,2)];
    H2O_HITRAN.data(:,2)=-log10(exp(-H2O_HITRAN.data(:,2)*10^-4));
    H2O_25T     =   [H2O_HITRAN.data(:,1),H2O_HITRAN.data(:,2)];
    
    [wave_start,CO2_25T_start_index]    =   ...
        min(abs(CO2_25T(:,1)-wavenumber_data(1)));
    [wave_end,CO2_25T_end_index]    =   ...
        min(abs(CO2_25T(:,1)-wavenumber_data(end)));

    [wave_start,H2O_25T_start_index]    =   ...
        min(abs(H2O_25T(:,1)-wavenumber_data(1)));
    [wave_end,H2O_25T_end_index]    =   ...
        min(abs(H2O_25T(:,1)-wavenumber_data(end)));
   
% These and further if's are used to flip the compound data vectors, which
% tend to have the highest k in (1,1), and the lowest in (end) for PNNL
% data.
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
        lower_H2O_25T_index     =   H2O_25T_start_index;-1;
        upper_H2O_25T_index     =   H2O_25T_end_index;+1;
    end
    
    CO2_25T_r     =   interp1( ...
        CO2_25T(lower_CO2_25T_index:upper_CO2_25T_index,1), ...
        CO2_25T(lower_CO2_25T_index:upper_CO2_25T_index,2), ...
        wavenumber_data);
    
    H2O_25T_r     =   interp1( ...
        H2O_25T(lower_H2O_25T_index:upper_H2O_25T_index,1), ...
        H2O_25T(lower_H2O_25T_index:upper_H2O_25T_index,2), ...
        wavenumber_data);


database_compounds_all(:,1) = CO2_25T_r;
database_compounds_all(:,2) = H2O_25T_r;

%allocate memory for the concentrations

concentrations = zeros(samples_healthy+samples_asthma:2);

for i=1:samples_healthy+samples_asthma;
    
    func_all = @(C) difference_for_least_squares_all( ...
        healthy_asthma(:,i), database_compounds_all*interaction_length, C);

    %Initial concentration guess in ppm.
        C_initial_all (1,1) = 1000;
        C_initial_all (2,1) = 1000;

    lb = zeros(size(C_initial_all)); %lower bound of zero
    ub = Inf*ones(size(C_initial_all));
    concentrations(i,:) = lsqnonlin(func_all,C_initial_all,lb,ub, options);
end

figure
% subplot(2,1,1);
plot(wavenumber_data,healthy_asthma(:,1),...
    wavenumber_data,interaction_length*bsxfun(@times,database_compounds_all',concentrations(1,:)'))
    xlabel('Wavenumber (cm^-^1)')
    ylabel('Absorbance')
    title('Absorbance calibrated by wavenumber');
legend('Absorbance','CO_2','H_2O')
% subplot(2,1,2);
% plot(wavenumber_data,healthy_asthma(:,2),...
%     wavenumber_data,interaction_length*bsxfun(@times,database_compounds_all',concentrations(2,:)'))
% legend('old absorbance','CO2','H2O')

