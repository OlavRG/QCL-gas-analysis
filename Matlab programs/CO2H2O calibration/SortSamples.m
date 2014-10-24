% INPUT
%     Sample Control_with removed entries.xlsx
%     20140703_Data_Processed_Wavenumber_Calibration.mat
%           Above is a MAT file with all measurement and monitor data, and
%           resulting wavenumber calibrated absorbances.

% FUNCTION
%    This script gets the absorbances from the .mat file and saves them to
%    the appropriate group and sample cells according to the xlsx file.

% importdata from L:\IST\OP\scratch\Olav Grouwstra\Sample
% Control_edited.xlsx as cell array 'sample_control', which tells what
% cells in Processed_Data_Wavenumber_Calibration belong to which group and
% sample. Some samples are left out on the advice of the people in
% Rotterdam taking breath samples, because the patients were on
% anti-biotics. The entries for these are deleted from 'Sample
% Control_edited.xlsx', which is the reason it's 'edited'.

addpath(genpath('L:\IST\OP\scratch\Olav Grouwstra\'));
cd('L:\IST\OP\scratch\Olav Grouwstra\');

[sam_num,sam_txt,sample_control] = xlsread(['L:\IST\OP\scratch\'...
    'Olav Grouwstra\Measurements\CO2H2O calibrated data\data of 26-8-2014\Sample Control_Checked_20140826_with Removed entries.xlsx']); 
load(['L:\IST\OP\scratch\Olav Grouwstra\Measurements\'...
    'Wavenumber calibrated data 3july2014\'...
    '20140703_Data_Processed_Wavenumber_Calibration.mat'])

G0_patients = (142-2)/4;
G1_patients = (80-2)/2;
G2_patients = (34-2)/2;
Absorbance_length = length(Processed_Data_Wavenumber_Calibration{1,1}.Absorbance);

G0S1=zeros(Absorbance_length,G0_patients);
G0S2=zeros(Absorbance_length,G0_patients);
G0S3=zeros(Absorbance_length,G0_patients);
G0S4=zeros(Absorbance_length,G0_patients);
G1S1=zeros(Absorbance_length,G1_patients);
G1S2=zeros(Absorbance_length,G1_patients);
G2S1=zeros(Absorbance_length,G2_patients);
G2S2=zeros(Absorbance_length,G2_patients);

% Separate m,n counters insure that empty cells in 'sample_control' won't be left
% as columns of zeros in the G#S# matrices.
m=1;n=1;
for k=1:G0_patients
    TUDsample1 = sample_control{-1+4*k,3};
    G0S1(:,k)=Processed_Data_Wavenumber_Calibration{1,TUDsample1}.Absorbance;

    TUDsample2 = sample_control{4*k,3};
    G0S2(:,k)=Processed_Data_Wavenumber_Calibration{1,TUDsample2}.Absorbance;
    
    TUDsample3 = sample_control{1+4*k,3};
    if isnan(TUDsample3)
    else
    G0S3(:,m)=Processed_Data_Wavenumber_Calibration{1,TUDsample3}.Absorbance;
    m=m+1;    
    end
    
    TUDsample4 = sample_control{2+4*k,3};
    if isnan(TUDsample4)
    else
    G0S4(:,n)=Processed_Data_Wavenumber_Calibration{1,TUDsample4}.Absorbance;
    n=n+1;    
    end
end
G0S3=G0S3(:,1:m-1);
G0S4=G0S4(:,1:n-1);

m=1;n=1;
for k=1:G1_patients
    TUDsample1 = sample_control{1+2*k,9};
    if isnan(TUDsample1)
    else
    G1S1(:,m)=Processed_Data_Wavenumber_Calibration{1,TUDsample1}.Absorbance;
    m=m+1;    
    end
    
    TUDsample2 = sample_control{2+2*k,9};
    if isnan(TUDsample2)
    else
    G1S2(:,n)=Processed_Data_Wavenumber_Calibration{1,TUDsample2}.Absorbance;
    n=n+1;    
    end
end
G1S1=G1S1(:,1:m-1);
G1S2=G1S2(:,1:n-1);

m=1;n=1;
for k=1:G2_patients
    TUDsample1 = sample_control{1+2*k,15};
    if isnan(TUDsample1)
    else
    G2S1(:,m)=Processed_Data_Wavenumber_Calibration{1,TUDsample1}.Absorbance;
    m=m+1;
    end
    
    TUDsample2 = sample_control{2+2*k,15};
    if isnan(TUDsample2)
    else
    G2S2(:,n)=Processed_Data_Wavenumber_Calibration{1,TUDsample2}.Absorbance;
    n=n+1;
    end
end
G2S1=G2S1(:,1:m-1);
G2S2=G2S2(:,1:n-1);

Wavenumber=Processed_Data_Wavenumber_Calibration{1,1}.Wavenumber;
save('Absorbances','Wavenumber','G0S1','G0S2','G0S3','G0S4','G1S1','G1S2','G2S1','G2S2')