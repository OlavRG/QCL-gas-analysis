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
function [ output1, output2, output3 ] = Acetone_and_H2O_Determination(input1, input2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

data_wavenumber = input1;
data_absorbance = input2;
interaction_length = 54.36; %Interaction length between laser and sample in meters

%% Load database files

if exist('Database_H2O_and_Acetone.mat','file') == 2
    %load the calibration file
    database_info = load('Database_H2O_and_Acetone.mat');
    database_wavenumber = database_info.Database_H2O_and_Acetone.wavenumber;
    database_H2O = database_info.Database_H2O_and_Acetone.H2O_base10_1ppm * interaction_length; % H2O database intensity for 1 ppm and the defined interaction length
    database_acetone = database_info.Database_H2O_and_Acetone.Acetone_base10_1ppm * interaction_length; % acetone database intensity for 1 ppm and the defined interaction length
else
    database_wavenumber = 0;
end

%% Part 1: Define function for least square minimization.

options = optimset('Display','iter','TolFun',1e-15); % 

func = @(C) difference_for_least_squares( data_absorbance, database_H2O, database_acetone, C );

C_initial = [1e3 1]; %Initial concentration guess. 1000 ppm for H2O and 1 ppm for acetone.  

%% Part 2: Optimize the concentrations. Use the least squared algorithm

[C_optimized,resnorm,residual] = lsqnonlin(func,C_initial,[],[], options);


%% Output
% 
 output1 = C_optimized;
 output2 = resnorm;
 output3 = residual;


end

