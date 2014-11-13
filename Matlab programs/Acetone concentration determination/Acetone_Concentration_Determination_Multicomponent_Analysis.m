function [ output1, output2, output3 ] = Acetone_Concentration_Determination_Multicomponent_Analysis( input1 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Dependencies
	% (calibrated) absorbance
	% Database_H2O_and_Acetone.mat
	% Calibration_H2O_v3
	% Acetone_and_H2O_Determination
	% Curve_Fitting_Removal_H2O
	% Multiple_Component_Determination
		% Database_Compounds_in_Breath_1150_1255.mat
		% difference_for_least_squares_all
        
% TODO
%     make script that gets input from 'data no zeros' text files and runs it.
%%

data = input1;
wavenumber = data.Wavenumber;
absorbance = data.Absorbance;

interaction_length = 54.36; 
%%

if exist('Database_H2O_and_Acetone.mat','file') == 2
    %load the calibration file
    database_info = load('Database_H2O_and_Acetone.mat');
    database_wavenumber = database_info.Database_H2O_and_Acetone.wavenumber;
    database_H2O = database_info.Database_H2O_and_Acetone.H2O_base10_1ppm * interaction_length; % H2O database intensity for 1 ppm and the defined interaction length
    database_acetone = database_info.Database_H2O_and_Acetone.Acetone_base10_1ppm * interaction_length; % acetone database intensity for 1 ppm and the defined interaction length
else
    database_wavenumber = 0;
end

%%

[ wavenumber_calibrated, absorbance_calibrated ] = Calibration_H2O_v3(wavenumber, absorbance);

[ concentrations_H2O_Acetone, resnorm_H2O_Acetone, residual_H2O_Acetone ] = Acetone_and_H2O_Determination(wavenumber_calibrated', absorbance_calibrated');

absorbance_No_Acetone = absorbance_calibrated - concentrations_H2O_Acetone(1,2)*database_acetone';

[ curve_fitting, absorbance_No_H2O ] = Curve_Fitting_Removal_H2O(wavenumber_calibrated', absorbance_No_Acetone');

absorbance_Acetone_back_No_H2O = absorbance_No_H2O + concentrations_H2O_Acetone(1,2)*database_acetone;

% plot(wavenumber_calibrated,absorbance_Acetone_back_No_H2O);
% hold on;

[  concentration_Acetone, resnorm_Acetone, residual_Acetone ] = Multiple_Component_Determination(absorbance_Acetone_back_No_H2O);


%%
 output1 = concentration_Acetone;
 output2 = resnorm_Acetone;
 output3 = residual_Acetone;


end

