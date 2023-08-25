function swr_cell_DFF = rip_cell_an(varargin)

% Written by Christoffer Berge

% Code that allows user to extract average signal from single ROIs during different 
% SWR types for multiple sessions and store data in a struct. For the first 
% session to be analyzed user only inputs "sData". For every following session 
% add "bulk_DFF" as a second argument and the code will update the struct with 
% the data from that sessions without overwriting the previous data. 

sData = varargin{1,1};
if length(varargin) > 1
    swr_cell_DFF = varargin{1,2};
end

%% Locate file name in folder (for saving file)
folder_name = dir;
folder_cell = strfind({folder_name.name}, '.mat');
for i = 1:length(folder_cell)
    if folder_cell{1,i} > 0
    idx = i;
    end
end

if length(idx) > 1
    error('Idx is > 1, could not locate sData file name')
end

fileName = folder_name(idx).name;

%% Run SWR bulk DF/F analysis
[rip_sig] = temp_code(sData);

names = {'_awake' ,'_awakeZ', '_single', '_singleZ', '_coupled', '_coupledZ'};
fileName = matlab.lang.makeValidName(fileName);
fileNames = strcat(fileName, names);
for i = 1:length(names)
    swr_cell_DFF.(fileNames{i}) = rip_sig{1,i};
end
clearvars -except sData swr_cell_DFF
