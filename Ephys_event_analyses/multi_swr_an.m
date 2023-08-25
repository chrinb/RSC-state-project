function [swr_mean_signal, opts] = multi_swr_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Code that allows user to get the mean DF/F during SWRs for multiple
% sessions and store data in a struct. For the first session to be analyzed
% user only inputs "sData". For every following session add the output (e.g.,
% "swr_mean_signal" as a second argument and "opts" as a third argument and 
% data from those sessions will be added to the struct using the same
% options. 

sData = varargin{1,1};

if length(varargin) > 1
    swr_mean_signal = varargin{1,2};
    opts            = varargin{1,3};
end

try isstruct(opts)

catch
    prompt = sprintf('Define initial settings to use  in analysis: ');
    
    settings    = input(prompt,'s'); %
    settings = split( string(settings));
    
    opts = get_multi_event_options(settings);
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
if opts.beh_state == 1
    rip_sig = swr_awake_an(sData, opts);
elseif opts.beh_state == 2
    rip_sig = swr_sleep_an(sData, opts);
end
% for sleep recordings
if opts.beh_state == 2
    names = {'_awake' ,'_awakeZ', '_single', '_singleZ', '_coupled', '_coupledZ'};
    fileName = matlab.lang.makeValidName(fileName);
    fileNames = strcat(fileName, names);
    for i = 1:length(names)
        swr_mean_signal.(fileNames{i}) = rip_sig{1,i};
    end
% for awake recordings
elseif opts.beh_state == 1
    names = {'_awake' ,'_awakeZ'};
    fileName = matlab.lang.makeValidName(fileName);
    fileNames = strcat(fileName, names);
    for i = 1:length(names)
        swr_mean_signal.(fileNames{i}) = rip_sig{1,i};
    end
end

% clearvars -except sData swr_mean_signal
