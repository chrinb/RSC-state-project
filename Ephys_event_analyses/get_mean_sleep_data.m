function [bouts_per_h, duration, total] = get_mean_sleep_data(varargin)

% Calculate mean 

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

% Index data
nrem_boutsH_idx   = (1:9:length(mean_subset_data));
rem_boutsH_idx    = (2:9:length(mean_subset_data));
is_boutsH_idx     = (3:9:length(mean_subset_data));
bout_nrem_dur_idx = (4:9:length(mean_subset_data));
bout_rem_dur_idx  = (5:9:length(mean_subset_data));
bout_is_dur_idx   = (6:9:length(mean_subset_data));
total_nrem_idx    = (7:9:length(mean_subset_data));
total_rem_idx     = (8:9:length(mean_subset_data));
total_is_idx      = (9:9:length(mean_subset_data));


boutsH_nrem_data   = mean_subset_data(nrem_boutsH_idx,:);
boutsH_rem_data     = mean_subset_data(rem_boutsH_idx,:);
boutsH_is_data     = mean_subset_data(is_boutsH_idx,:);

bout_nrem_dur_data     = mean_subset_data(bout_nrem_dur_idx,:);
bout_rem_dur_data   = mean_subset_data(bout_rem_dur_idx,:);
bout_is_dur_data   = mean_subset_data(bout_is_dur_idx,:);

total_nrem_data   = mean_subset_data(total_nrem_idx,:);
total_rem_data   = mean_subset_data(total_rem_idx,:);
total_is_data   = mean_subset_data(total_is_idx,:);

% Convert from cell to array
boutsH_nrem_array   = cell2mat(boutsH_nrem_data);
boutsH_rem_array     = cell2mat(boutsH_rem_data);
boutsH_is_array    = cell2mat(boutsH_is_data);

bout_nrem_dur_array     = cell2mat(bout_nrem_dur_data);
bout_rem_dur_array   = cell2mat(bout_rem_dur_data);
bout_is_dur_array   = cell2mat(bout_is_dur_data);

total_nrem_array     = cell2mat(total_nrem_data);
total_rem_array     = cell2mat(total_rem_data);
total_is_array     = cell2mat(total_is_data);


bouts_per_h = {boutsH_nrem_array, boutsH_rem_array, boutsH_is_array};
duration    = {bout_nrem_dur_array, bout_rem_dur_array,bout_is_dur_array};
total       = {total_nrem_array,total_rem_array,total_is_array};
