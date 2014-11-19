function [ F ] = difference_for_least_squares_all ...
    ( data_absorbance, database_compounds_all, C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%[compounds_points, num_compounds] = size(database_compounds_all);
model = database_compounds_all*C;
F = data_absorbance - model;

end

