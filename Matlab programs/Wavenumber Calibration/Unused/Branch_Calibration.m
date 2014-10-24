function [ output1 ] = Branch_Calibration(input1, input2, input3 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

wavenumber_database = input1(:,1);
data_wavenumber = input2(:,1);
data_absorbance = input2(:,2);
zone = input3;


%% Find the closest peak in the data to the reference peaks in the database

for i = 1:length(wavenumber_database)
   
    index(i) =  find(data_wavenumber >= wavenumber_database(i,1),1,'first');
    start_search_index = index(i) - 5;
    end_search_index = index(i) + 15;
    search_zone = data_absorbance(start_search_index:end_search_index);
    max_search= max(search_zone);
    max_index(i) = find(data_absorbance == max_search,1,'first');
    
end


%% Linear fit

fitting_ref= wavenumber_database;
fitting_meas= data_wavenumber(max_index);
cuadfit = fit(fitting_meas,fitting_ref,'poly1');

% plot (fitting_meas,fitting_ref,'r*');
% hold on

if zone == 0
    data_wavenumber_zone = data_wavenumber(max_index(1):max_index(end));
    data_absorbance_zone = data_absorbance(max_index(1):max_index(end));
elseif zone == 1
    data_wavenumber_zone = data_wavenumber(max_index(3)+1:max_index(end-2)-1);
    data_absorbance_zone = data_absorbance(max_index(3)+1:max_index(end-2)-1);
elseif zone == 2
    data_wavenumber_zone = data_wavenumber(max_index(1):end);
    data_absorbance_zone = data_absorbance(max_index(1):end);
end

calibrated_wavenumber = feval(cuadfit,data_wavenumber_zone);

calibrated_data(:,1) = calibrated_wavenumber;
calibrated_data(:,2) = data_absorbance_zone;

% plot(calibrated_data(:,1),calibrated_data(:,2),'b');
% hold on


output1 = calibrated_data;

end

