function [ F ] = difference_for_least_squares( data_absorbance, database_H2O, database_acetone, C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
k = 1:length(database_H2O);
F = data_absorbance(k) - (C(1)*database_H2O(k) + C(2)*database_acetone(k));
end

