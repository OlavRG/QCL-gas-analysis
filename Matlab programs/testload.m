clear 'filestring' 'fid' 'tmp' 

tic
for k=1:5
filestring='L:\IST\OP\scratch\Olav Grouwstra\Compounds\Carbon_monoxide\CO_25T.TXT';

% 0.39
% tmp=load(filestring);

% 0.27
% tmp=importdata(filestring);

% 0.192
% fid = fopen(filestring);
% tmp = fscanf(fid, '%f %f');
% fclose(fid);

% 0.15
% fid = fileread(filestring);
% tmp = textscan(fid, '%f %f');
 
% 0.13
fid = fopen(filestring);
tmp = textscan(fid, '%f %f');
fclose(fid);

time=toc;
disp([num2str(time) ' ' num2str(k) ' out of 5']) 
end