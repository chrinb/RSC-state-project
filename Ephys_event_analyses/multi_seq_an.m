function mean_event_rate = multi_seq_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Code that allows user to extract nr of delta-spindle, SWR-delta, and
% SWR-delta-spindle events pr NREM minute across multiple sessions. For the
% first session only input "sData". For subsequent sessions, include
% "mean_event_rate" output as a second argument to the function to prevent
% overwriting the values from the previous analyzed session. 

sData = varargin{1,1};

if length(varargin) > 1
    mean_event_rate = varargin{1,2};
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

oscillation_seq = oscillation_coupling(sData);
names = {'_delta_spindle', '_SWR_delta', '_SWR_delta_spindle'};
fileName = matlab.lang.makeValidName(fileName);
fileNames = strcat(fileName, names);
for i = 1:length(names)
    mean_event_rate.(fileNames{i}) = oscillation_seq{1,i};
end