    15-10-2014 by Olav Grouwstra
	
	70 healthy samples, 70 asthmatic samples
	
	max_p           =   0.05;
    min_p_region    =   10;
    min_intensity   =   10^(-4);
    interaction_length = 54.36; % in meters

    options = optimset('Display','off','TolFun',1e-15); 
    region_pres=1;
