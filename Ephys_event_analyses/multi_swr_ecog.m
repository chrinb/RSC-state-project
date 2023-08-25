function swr_mean_signal = multi_swr_ecog(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Code that allows user to get the mean SWR-aligned ECoG activity for multiple
% sessions and store data in a struct. For the first session to be analyzed
% user only inputs "sData". For every following session add the output (e.g.,
% "swr_mean_signal" as a second argument and the code will update the struct 
% with the data from that sessions without overwriting the previous data. 

sData = varargin{1,1};
if length(varargin) > 1
    swr_mean_signal = varargin{1,2};
end

% prompt = sprintf('Sleep (1) or awake(2) recordings? ');
% awake_or_sleep = input(prompt);

% By default set to sleep. Change for awake recordings.
awake_or_sleep = 1;


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
if awake_or_sleep == 1
    rip_sig = swr_aligned_ecog(sData, 1 , [] , 4 , 1,[]);
    names = {'_unclassified' ,'_unclassifiedZ', '_single', '_singleZ', '_coupled', '_coupledZ'};
    fileName = matlab.lang.makeValidName(fileName);
    fileNames = strcat(fileName, names);
    for i = 1:length(names)
        swr_mean_signal.(fileNames{i}) = rip_sig{1,i};
    end
elseif awake_or_sleep == 2
    rip_sig = swr_aligned_ecog_awake(sData);
    names = {'_awake' ,'_awakeZ'};
    fileName = matlab.lang.makeValidName(fileName);
    fileNames = strcat(fileName, names);
    for i = 1:length(names)
        swr_mean_signal.(fileNames{i}) = rip_sig{1,i};
    end
end


clearvars -except sData swr_mean_signal
