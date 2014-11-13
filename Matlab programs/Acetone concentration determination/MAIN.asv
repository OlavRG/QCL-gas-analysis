

%   Input
data_path   =   ['D:\Workspace\Breath Analysis\Measurements\data of 26-8-2014\data no zeros\'];
list        =   ['G0S1';'G0S2';'G0S3';'G0S4';'G1S1';'G1S2';'G2S1';'G2S2'];
G0S1        =   importdata([data_path list(1,:) '.TXT']); G0S1=G0S1.data;
G0S2        =   importdata([data_path list(2,:) '.TXT']); G0S2=G0S2.data;
G0S3        =   importdata([data_path list(3,:) '.TXT']); G0S3=G0S3.data;
G0S4        =   importdata([data_path list(4,:) '.TXT']); G0S4=G0S4.data;
G1S1        =   importdata([data_path list(5,:) '.TXT']); G1S1=G1S1.data;
G1S2        =   importdata([data_path list(6,:) '.TXT']); G1S2=G1S2.data;
G2S1        =   importdata([data_path list(7,:) '.TXT']); G2S1=G2S1.data;
G2S2        =   importdata([data_path list(8,:) '.TXT']); G2S2=G2S2.data;
load('D:\Workspace\Breath Analysis\Measurements\Wavenumber.txt')


healthy     =   [G0S1,G0S2];
asthma      =   [G1S1,G1S2];

healthy_concentration_Acetone = zeros(1,70);
healthy_resnorm_Acetone = zeros(1,70);
healthy_residual_Acetone = zeros(2101,70);

matlabpool local 2
parfor k=1:70;
    data = [];
    data.Wavenumber = Wavenumber; 
    data.Absorbance = healthy(:,k);
    [ healthy_concentration_Acetone(1,k), healthy_resnorm_Acetone(1,k), ...
        healthy_residual_Acetone(:,k) ] = ...
        Acetone_Concentration_Determination_Multicomponent_Analysis( data );
    k
end

asthma_concentration_Acetone = zeros(1,70);
asthma_resnorm_Acetone = zeros(1,70);
asthma_residual_Acetone = zeros(2101,70);

parfor k=1:70;
    data = [];
    data.Wavenumber = Wavenumber;
    data.Absorbance = asthma(:,k);
    [ asthma_concentration_Acetone(1,k), asthma_resnorm_Acetone(1,k), ...
        asthma_residual_Acetone(:,k) ] = ...
        Acetone_Concentration_Determination_Multicomponent_Analysis( data );
    k
end
matlabpool close