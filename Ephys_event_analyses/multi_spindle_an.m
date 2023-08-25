function spindle_mean_signal = multi_spindle_an(varargin)

% Written by Christoffer Berge

% Code that allows user to get the mean DF/F during sleep spindles for multiple
% sessions and store data in a struct. For the first session to be analyzed
% user only inputs "sData". For every following session add "bulk_DFF" as a
% second argument and the code will update the struct with the data from
% that sessions without overwriting the previous data. 

sData = varargin{1,1};
if length(varargin) > 1
    spindle_mean_signal = varargin{1,2};
end

% prompt = sprintf('Bulk (1) or single ROI (2) data? ');
% select = input(prompt);

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

%% Run mean spindle activity analysis
% if select == 1
%     [~,~,spindle_sig] = spinSigAn2(sData);
% elseif select == 2
%     spindle_sig = spinSigAn(sData);
% end
spindle_sig = rand_nrem_time_comp(sData);

names = {'_start' ,'_startZ'}; % raw and z-scored data
fileName = matlab.lang.makeValidName(fileName);
fileNames = strcat(fileName, names);
for i = 1:length(names)
    spindle_mean_signal.(fileNames{i}) = spindle_sig{1,i};
end
clearvars -except sData spindle_mean_signal
