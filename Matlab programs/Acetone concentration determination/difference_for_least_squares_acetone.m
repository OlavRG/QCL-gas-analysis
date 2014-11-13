function [ F ] = difference_for_least_squares_acetone( data_absorbance, database_acetone, C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
k = 1:length(database_acetone);
lin = k;
F = data_absorbance(k) - (C(1)*database_acetone(k));
end

