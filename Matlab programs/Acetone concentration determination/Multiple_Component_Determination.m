%% ACETONE_AND_H2O_DETERMINATION
% C
% Syntax:
%
%       [
% 
% Input:
%         
%       
%               
% Output:
%   
%       
%
% Example(s):
%        
%      
% Dependence: NONE
% 
% Author: Adonis Reyes Reyes
% Email: adonisatl@gmail.com
% 
% Project: IR Spectroscopy
% 
%
%Version: v1.0 (18-02-2014)
%
%% Change-log
%
% v1.0 (20-02-2014)  First operation version - Calibration by parts
 
%% Function declaration
function [ output1, output2, output3] = Multiple_Component_Determination(input1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

data_absorbance = input1;
interaction_length = 54.36; %Interaction length between laser and sample in meters

%% Load database files

if exist('Database_Compounds_in_Breath_1150_1255.mat','file') == 2
    %load the calibration file
    database_info = load('Database_Compounds_in_Breath_1150_1255.mat');
    database_wavenumber = database_info.Database_Compounds_in_Breath_1150_1255.wavenumber;
    
    database_compounds_all(:,1) = database_info.Database_Compounds_in_Breath_1150_1255.Acetone_base10_1ppm;
%     database_compounds_all(:,2) = database_info.Database_Compounds_in_Breath_1150_1255.CO2_base10_1ppm;
%     database_compounds_all(:,3) = database_info.Database_Compounds_in_Breath_1150_1255.Acetic_Acid_base10_1ppm;
%     database_compounds_all(:,4) = database_info.Database_Compounds_in_Breath_1150_1255.Ammonia_base10_1ppm;
%     database_compounds_all(:,5) = database_info.Database_Compounds_in_Breath_1150_1255.Methane_base10_1ppm;
%    database_compounds_all(:,6) = database_info.Database_Compounds_in_Breath_1150_1255.H2O_base10_1ppm;
%    database_compounds_all(:,7) = database_info.Database_Compounds_in_Breath_1150_1255.Ethanol_base10_1ppm;
%    database_compounds_all(:,8) = database_info.Database_Compounds_in_Breath_1150_1255.Formaldehyde_base10_1ppm;
%    database_compounds_all(:,9) = database_info.Database_Compounds_in_Breath_1150_1255.Pentane_base10_1ppm;
%    database_compounds_all(:,10) = database_info.Database_Compounds_in_Breath_1150_1255.Benzene_base10_1ppm;
    
else
    database_wavenumber = 0;
end

% a=find(data_wavenumber == 1150, 1,'first');
% b=find(data_wavenumber == 1255, 1, 'first');
% 
% absorbance_region = data_absorbance(a:b);


%% Part 1: Define function for least square minimization.

options = optimset('Display','iter','TolFun',1e-15); % 

    func_all = @(C) difference_for_least_squares_all( ...
        data_absorbance, database_compounds_all*interaction_length, C);
    
    %Initial concentration guess in ppm.
        C_initial_all (1,1) = 1;
%         C_initial_all (2,1) = 1000;
%         C_initial_all (3:5,1) = 0.1;
%         C_initial_all (6,1) = 1000;
        

    lb = zeros(size(C_initial_all)); %lower bound of zero
    ub = Inf*ones(size(C_initial_all));
    [C_optimized,resnorm,residual] = lsqnonlin(func_all,C_initial_all,lb,ub, options);
        
%func = @(C) difference_for_least_squares( data_absorbance, database_H2O, database_acetone, C );

%C_initial = [1e3 1]; %Initial concentration guess. 1000 ppm for H2O and 1 ppm for acetone.  

%% Part 2: Optimize the concentrations. Use the least squared algorithm

%[C_optimized,resnorm,residual] = lsqnonlin(func,C_initial,[],[], options);


%% Output
% 
 output1 = C_optimized;
 output2 = resnorm;
 output3 = residual;


end

