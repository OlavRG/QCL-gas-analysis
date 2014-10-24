
% FUNCTION
%     Removes all rows from G0S1 not in the ind variable, and writes the
%     new rows to new files. The ind variable has only those rows with low
%     p-values (<0.005).

% INPUT
%   ind variable from HealhtyAsthma_analysis.m
%   G0S1.txt
tic;

dataPath = ['L:\IST\OP\scratch\Olav Grouwstra\Measurements\'...
    'CO2H2O calibrated data\data of 26-8-2014\data no zeros\'];
files=dir(dataPath);

for k=3:10
    filename = [dataPath files(k,1).name];
    delimiter = '\t';
    startRow = 2;
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    GXSX = [dataArray{1:end-1}];
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;

    GXSX = GXSX(ind,:);

    dlmwrite([dataPath 'lim' files(k,1).name],GXSX,...
    '-append','delimiter','\t','newline','pc','precision','%2.16f')
    if k==8 halt 
    end
    toc
end