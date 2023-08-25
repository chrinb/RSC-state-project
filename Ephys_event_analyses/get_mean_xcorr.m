function data_output = get_mean_xcorr(varargin)

% Calculate mean cross-correlation for different events (SWRs, sleep spindles,
% slow waves). First input argument
% has to be a struct containing the mean data cross correlations from 
% several sessions. If sessions from multiple
% mice are loaded in struct, specify the mouse name ("m6120"). 

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
swr_spindle_idx   = (1:5:length(mean_subset_data));
swr_delta_idx     = (2:5:length(mean_subset_data));
delta_spindle_idx = (3:5:length(mean_subset_data));
delta_swr_idx     = (4:5:length(mean_subset_data));
spindle_swr_idx   = (5:5:length(mean_subset_data));

swr_spindle_data   = mean_subset_data(swr_spindle_idx,:);
swr_delta_data     = mean_subset_data(swr_delta_idx,:);
delta_spindle_data = mean_subset_data(delta_spindle_idx,:);
delta_swr_data     = mean_subset_data(delta_swr_idx,:);
spindle_swr_data   = mean_subset_data(spindle_swr_idx,:);

% Convert from cell to array
swr_spindle_array   = cell2mat(swr_spindle_data);
swr_delta_array     = cell2mat(swr_delta_data);
delta_spindle_array = cell2mat(delta_spindle_data);
delta_swr_array     = cell2mat(delta_swr_data);
spindle_swr_array   = cell2mat(spindle_swr_data);

data_output = {swr_spindle_array, swr_delta_array, delta_spindle_array,...
    delta_swr_array, spindle_swr_array};