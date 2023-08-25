function multi_sleep_feat = multi_sleep_feat(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Code that allows user to extract nr of delta-spindle, SWR-delta, and
% SWR-delta-spindle events pr NREM minute across multiple sessions. For the
% first session only input "sData". For subsequent sessions, include
% "mean_event_rate" output as a second argument to the function to prevent
% overwriting the values from the previous analyzed session. 

sData = varargin{1,1};

if length(varargin) > 1
    multi_sleep_feat = varargin{1,2};
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

mean_sleep_feat = sleepfeatures(sData);
names = {'_boutsNREM_perH', '_boutsREM_perH', '_boutsIS_perH', '_boutNREM_dur',...
     '_boutREM_dur', '_boutIS_dur', '_totalNREMmin', '_totalREMmin', '_totalISmin'};
fileName = matlab.lang.makeValidName(fileName);
fileNames = strcat(fileName, names);
for i = 1:length(names)
    multi_sleep_feat.(fileNames{i}) = mean_sleep_feat{1,i};
end
