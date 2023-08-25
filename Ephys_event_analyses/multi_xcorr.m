function mean_multi_xcorr = multi_xcorr(varargin)


% Written by Christoffer Berge | Vervaeke Lab

% Code that allows user to extract mean cross-correlations across multiple 
% sessions using the "event_cross_corr" function. For the
% first session only input "sData". For subsequent sessions, include
% "mean_multi_xcorr" output as a second argument to the function to prevent
% overwriting the values from the previous analyzed session. 

sData = varargin{1,1};

if length(varargin) > 1
    mean_multi_xcorr = varargin{1,2};
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

%% Get oscillation sequence data

mean_xcorr = event_cross_corr(sData);
names = {'_swr_spindle', '_swr_delta', '_delta_swr', '_delta_spindle', '_spindle_swr'};
fileName = matlab.lang.makeValidName(fileName);
fileNames = strcat(fileName, names);
for i = 1:length(names)
    mean_multi_xcorr.(fileNames{i}) = mean_xcorr{1,i};
end