function [ output1 ] = Calibration(input1)


%% Load calibration file
if exist('WavenumberCalibration.mat','file') == 2
    %load the calibration file
    load('WavenumberCalibration.mat')
else
    wavenumber_database = 0;
end



data(:,1) = input1(:,1);
data(:,2) = input1(:,2);


%% Copy first part of the spectrum

end_zone0 = find (data(:,1) == Wavenumber_Ref.zone1(1,1));
calibrated_data(:,1) = data(1:end_zone0-1,1);
calibrated_data(:,2) = data(1:end_zone0-1,2);

%% Calibrate each branch

calibrated_data = [calibrated_data ; Branch_Calibration(Wavenumber_Ref.zone1(:,1), data, 0 )];
calibrated_data = [calibrated_data ; Branch_Calibration(Wavenumber_Ref.zone2(:,1), data, 1 )];
calibrated_data = [calibrated_data ; Branch_Calibration(Wavenumber_Ref.zone3(:,1), data, 0 )];
calibrated_data = [calibrated_data ; Branch_Calibration(Wavenumber_Ref.zone4(:,1), data, 1 )];
calibrated_data = [calibrated_data ; Branch_Calibration(Wavenumber_Ref.zone5(:,1), data, 2 )];

wavenumber_interpolation = [832 : 0.05 : 1263]';
absorbance_interpolation = interp1(calibrated_data(:,1),calibrated_data(:,2),wavenumber_interpolation);

% interpolated_data(:,1) = wavenumber_interpolation;
interpolated_data(:,1) = absorbance_interpolation;

%% Function Part3: Remove NaN values

interpolated_data= interpolated_data(~any(isnan(interpolated_data), 2), :);

%  figure
% plot(calibrated_data(:,1),calibrated_data(:,2),'b');
% hold on
% plot(wavenumber_interpolation,vq,'r');

output1 = interpolated_data;

end

