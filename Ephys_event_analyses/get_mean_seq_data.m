function [nrem_delta_spindle,nrem_swr_delta,nrem_swr_delta_spindle] = ...
    get_mean_seq_data(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Calculate mean nr of delta-spindle, SWR-delta, and SWR-delta-spindle 
% events across multiple sessions. Input has to be a struct containing the
% NREM event rate for the three events organized as 1 = delta-spindle session 1
%, 2 = SWR-delta session 1, 3 = SWR-delta-spindle session 1, 4 = 
% delta-spindle session 2etc. Use in conjuction with function
% "multi_seq_an". 

data_struct = varargin{1,1};
fieldName = varargin{1,2}; 
if size(varargin,2) == 3
    fieldName2 = varargin{1,3};
end


% Convert struct to cell for indexing
data_cell = struct2cell(data_struct);

% Get field names in struct
f = fieldnames(data_struct(1));

% Initialize data matrix
mean_subset_data = [];

% Loop over field names 
for i=1:length(f)
    % if two field names are specified, check that each field name contains
    % both strings
    if size(varargin,2) == 3
        split_field = split(f{i}, "_");
        if contains( split_field(1), fieldName) && contains( split_field(5), fieldName2)
            mean_subset_data = vertcat(mean_subset_data, data_cell(i) );
        end

    elseif contains( f{i}, strtrim(fieldName),'IgnoreCase',true)
        mean_subset_data = vertcat(mean_subset_data, data_cell(i) );
    end
end

% Index data and z-scored data separately
delta_spindle_idx     = (1:3:length(mean_subset_data));
swr_delta_idx         = (2:3:length(mean_subset_data));
swr_delta_spindle_idx = (3:3:length(mean_subset_data));

nrem_delta_spindle_cell     = mean_subset_data(delta_spindle_idx,:);
nrem_swr_delta_cell         = mean_subset_data(swr_delta_idx,:);
nrem_swr_delta_spindle_cell = mean_subset_data(swr_delta_spindle_idx,:);
% Convert from cell to array
nrem_delta_spindle     = cell2mat(nrem_delta_spindle_cell);
nrem_swr_delta         = cell2mat(nrem_swr_delta_cell);
nrem_swr_delta_spindle = cell2mat(nrem_swr_delta_spindle_cell);
