function multi_swr_sleep = multi_swr_feat_sleep(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Code that allows user to extract various NREM SWR features across multiple 
% sessions. For the first session only input "sData". For subsequent sessions,
% include as a second input argument the output of the first time the function
% was ran. E.g., input "multi_swr_sleep" for the second and all following iterations
% of the function 

sData = varargin{1,1};

if length(varargin) > 1
    multi_swr_sleep = varargin{1,2};
end

%% Locate session data file name in folder 
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

%% Get sleep session SWR data

swr_nrem_feat = swr_features_sleep(sData);

variable_names_string_array = string(swr_nrem_feat(1,:));
variable_names              = matlab.lang.makeValidName(variable_names_string_array);
% names         = {'_abs_swr_nr_in_rec',         '_proportion_clust_nrem_swr', ...
%                  '_proportion_clust_swr',      '_nrem_swr_dur_in_sec', ...
%                  '_nrem_swr_rate_per_min',     '_zScoreAmp',...
%                  '_nr_spindle_coupled_swrs',   '_proportion_spindle_coupled_swrs', ...
%                  '_nr_swr_clusters_in_rec',    '_prop_swr_clusters_in_rec' , ...
%                  '_nr_swr_cluster_in_nrem',    '_prop_swr_cluster_in_nrem' ,...
%                  '_nr_swr_cluster_in_spindle', '_prop_swr_cluster_in_spindle'...
%                  '_ISI'};
session_name = matlab.lang.makeValidName(session_name);
fileNames    = strcat(session_name, '_', variable_names);
for i = 1:length(variable_names)
    multi_swr_sleep.(fileNames{i}) = swr_nrem_feat{2,i};
end
