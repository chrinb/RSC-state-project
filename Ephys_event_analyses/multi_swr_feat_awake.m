function multi_swr_awake = multi_swr_feat_awake(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Code that allows user to extract various SWR features across multiple 
% *awake sessions. For the first session only input "sData". For subsequent sessions,
% include as a second input argument the output of the first time the function
% was ran. E.g., input "multi_swr_awake" for the second and all following iterations
% of the function 

sData = varargin{1,1};

if length(varargin) > 1
    multi_swr_awake = varargin{1,2};
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

session_name = folder_name(idx).name;

%% Get session SWR data

swr_awake_feat = swr_features_awake(sData);

variable_names_string_array = string(swr_awake_feat(1,:));
variable_names              = matlab.lang.makeValidName(variable_names_string_array);
% names = {'_abs_swr_nr_in_rec', '_awake_swr_rate_per_min', '_proportion_clust_swr', ...
%     '_swr_dur_in_sec', '_zScoreAmp', '_nr_swr_clusters_in_rec', ...
%     '_prop_swr_clusters_in_rec', '_ISI'};
session_name = matlab.lang.makeValidName(session_name);
fileNames = strcat(session_name, '_',variable_names);
for i = 1:length(variable_names)
    multi_swr_awake.(fileNames{i}) = swr_awake_feat{2,i};
end
