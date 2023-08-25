function [data_all_mice, miceIDs] = load_multi_session_data(varargin)

% Written by Eivind Hennestad. Modified by Christoffer Berge


% Function to load multiple datafiles into struct. 

prompt = sprintf('Type mouse ID: '); % user types in mouse ID, e.g. "m6120"

miceIDs = input(prompt, 's');
miceIDs = split(miceIDs, ' ')';

% Find absolute path to data
path = cell(size(miceIDs,2), 1 );
% Loop over nr of mice
for i = 1:size(miceIDs,2)
    full_miceIDs = replaceBetween(miceIDs{1,i}, 'm', '6', 'ouse');
    path{i}      = strcat('D:\Data\', convertCharsToStrings(full_miceIDs), '\Results\');
end



% Load data
data_all_mice = struct();
if strcmp(varargin{1,1}{5}, 'avg')
    if strcmp(varargin{1, 1}{1, 2}, 'swr')
        txt1 = 'mean_swr_';
    elseif strcmp(varargin{1, 1}{1, 2}, 'spindle')
        txt1 = 'mean_spindle';
    elseif strcmp(varargin{1, 1}{1, 2}, 'swa')
        txt1 = 'mean_swa';
    end
elseif strcmp(varargin{1,1}{5}, 'mod')
    txt1 = 'swr_mod';
end

% Check if user has specified string to select a particular dataset in
% Results folder. The defaulet is _mean_swr_ ... .mat
if size(varargin{1, 1},2) > 4
    specify_dataset = varargin{1, 1}{6};
    data_str = [txt1, varargin{1,1}{4},'_',varargin{1,1}{1}, specify_dataset,'.mat'];
else
    data_str = [txt1, varargin{1,1}{4},'_',varargin{1,1}{1} '.mat'];
end
% Loop over nr of mice
for i = 1:numel(miceIDs)
    % Select absolute path to data
    data_path = [ path{i},  [miceIDs{1,i},data_str ] ];
    data_path = append(data_path(1), data_path(2) );
    % Load data into struct
    data_all_mice.(miceIDs{i}) = load(data_path );
end