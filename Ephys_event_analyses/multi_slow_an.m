function slow_mean_signal = multi_slow_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Code that allows user to get the mean DF/F during slow wave activity for multiple
% sessions and store data in a struct. For the first session to be analyzed
% user only inputs "sData". For every following session add the output (e.g.,
% "slow_mean_signal" as a second argument and the code will update the struct 
% with the data from that sessions without overwriting the previous data. 

sData = varargin{1,1};
if length(varargin) > 1
    slow_mean_signal = varargin{1,2};
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

%% Run slow wave analysis

slow_sig =  so_sig_an(sData);

% for sleep recordings
names = {'_SO' , '_SOz', '_delta', '_deltaZ'};
fileName = matlab.lang.makeValidName(fileName);
fileNames = strcat(fileName, names);
for i = 1:length(names)
    slow_mean_signal.(fileNames{i}) = slow_sig{1,i};
end

clearvars -except sData slow_mean_signal
