function [data_all_mice, miceIDs] = load_multi_session_state_activity(varargin)

% Written by Christoffer Berge | Vervake lab

% Load cross-correlation data for plot function

prompt = sprintf('Type mouse ID: '); % user types in mouse ID, e.g. "m6120"

miceIDs = input(prompt, 's');
miceIDs = split(miceIDs, ' ')';

% Find absolute path to data
path = cell(size(miceIDs,2), 1 );

% Loop over nr of mice
for i = 1:size(miceIDs,2)
    full_miceIDs = replaceBetween(miceIDs{1,i}, 'm', '6', 'ouse');
    path{i}      = strcat('E:\Data\', convertCharsToStrings(full_miceIDs), '\Results\');
end

% Load data
data_str = '_state_activity.mat';

% Loop over nr of mice
for i = 1:numel(miceIDs)
    % Select absolute path to data
    data_path = [ path{i} miceIDs(i) data_str ];
    data_path = append(data_path(1), data_path(2), data_path(3) );
    % Load data into struct
    data_all_mice.(miceIDs{i}) = load(data_path );
end