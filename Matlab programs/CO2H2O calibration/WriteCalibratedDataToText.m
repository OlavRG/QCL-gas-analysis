
% FUNCTION
%     This file converts adjustedH2OCO2_absorbance to G0S1.txt etc. files.
%     This can be used to 
% 

% INPUT
%     adjustedH2OCO2_absorbance
%     G0S1 etc. matrices
%         only used for size determination 
%     Sample Control_with removed entries.xlsx

%     alternative for adjustedH2OCO2_absorbance: G0S1 etc. matrices 
%       (Need modification of the dlmwrite to get the data from such
%       matrices.)


% OUTPUT
%     G0S1.txt etc. files with sample numbers (e.g. 023/02) in first row


function [filename]=WriteCalibratedDataToText(adjustedH2OCO2_absorbance,...
    sampleControlPath,saveTXTsTo,G0S1,G0S2,G0S3,G0S4,G1S1,G1S2,G2S1,G2S2)
tic
[sam_num,sam_txt,sample_control] = xlsread(sampleControlPath); 

list={'G0S1'
'G0S2'
'G0S3'
'G0S4'
'G1S1'
'G1S2'
'G2S1'
'G2S2'
};
[lengthlist,sampleslist]    =   size(list);

list2={G0S1
G0S2
G0S3
G0S4
G1S1
G1S2
G2S1
G2S2
};
samples_group   =   [4,4,4,4,2,2,2,2];
start_column    =   [2,2,2,2,8,8,14,14];
start_row       =   [3,4,5,6,3,4,3,4];

start_sample=1;
end_sample=0;
for k=1:lengthlist
    
    
    [rows, sample_amount] =   size(list2{k});
    end_sample              =   end_sample + sample_amount;
    filename                =   strcat(list(k),'.txt');
    sample_numbers          =   sample_control(...
        start_row(k):samples_group(k):start_row(k)+(sample_amount-1)*samples_group(k),...
        start_column(k));
    fid=fopen([saveTXTsTo filename{1}],'wt');
    [rows,cols]=size(sample_numbers);
      fprintf(fid,'%s\t',sample_numbers{1:end-1,1:cols});
      fprintf(fid,'%s\n',sample_numbers{end,1:cols});
    fclose(fid);

%     dlmwrite(filename{1},sample_numbers,'delimiter','\t','newline','pc','precision','%2.16f')
    dlmwrite([saveTXTsTo filename{1}],adjustedH2OCO2_absorbance(:,start_sample:end_sample),...
        '-append','delimiter','\t','newline','pc','precision','%2.16f')
    
    start_sample   =   start_sample + sample_amount;
toc
end
end