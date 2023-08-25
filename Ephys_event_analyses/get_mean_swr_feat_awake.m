function mean_feat_awake = get_mean_swr_feat_awake(varargin)

% Calculate mean 

multi_session_struct = varargin{1,1};
% fieldName = varargin{1,2}; 
% if size(varargin,2) == 3
%     fieldName2 = varargin{1,3};
% end


% Convert struct to cell for indexing
data_cell = struct2cell(multi_session_struct);
mean_subset_data = data_cell;

% Get field names in struct
field_names_cell = fieldnames(multi_session_struct(1));


% Convert field names from multi-session struct (which is stored as a cell array)
% to string array, split string array at mouse+session ID | variable names, find unique
% variable names and nr of unique variable names
all_fields_string_array  = string(field_names_cell);
string_array_split       = split(all_fields_string_array, '_mat_');
all_variable_names       = string_array_split(:,2);
unique_variable_names    = unique(all_variable_names, 'stable');
nr_of_unique_field_names = length(unique_variable_names);


% Create a new cell array where first row contains unique variable names
mean_feat_awake = cellstr(unique_variable_names');

% Loop over nr of unique variable names 
for i = 1:nr_of_unique_field_names

    % The data in the multi-session struct repeats every 
    % "nr-of-unique-variable-names"-th time. Use this nr to index repeating
    % variables. 
    data_idx             = i:nr_of_unique_field_names:length(mean_subset_data);
    data_cell            = mean_subset_data(data_idx,:);
    data_arr             = cell2mat(data_cell);
    mean_feat_awake{2,i} = data_arr;
end



% Initialize data matrix
% mean_subset_data = [];

% % Loop over field names 
% for i=1:length(field_names_cell)
%     % if two field names are specified, check that each field name contains
%     % both strings
%     if size(varargin,2) == 3
%         split_field = split(field_names_cell{i}, "_");
%         if contains( split_field(1), fieldName) && contains( split_field(5), fieldName2)
%             mean_subset_data = vertcat(mean_subset_data, data_cell(i) );
%         end
% 
%     elseif contains( field_names_cell{i}, strtrim(fieldName),'IgnoreCase',true)
%         mean_subset_data = vertcat(mean_subset_data, data_cell(i) );
%     end
% end
% 
% % Index data
% absSWRinRec_idx            = (1:9:length(mean_subset_data));
% awake_swr_rate_per_min_idx = (2:5:length(mean_subset_data));
% propClustSWR_idx           = (3:5:length(mean_subset_data));
% swr_dur_in_sec_idx         = (4:5:length(mean_subset_data));
% swrAmp_idx                 = (5:5:length(mean_subset_data));
% nr_swr_clusters_in_rec     = (6:6:length(mean_subset_data));
% prop_swr_clusters_in_rec   = (7:8:length(mean_subset_data));
% ISI                        = (8:9:length(mean_subset_data));
% 
% absSWRinRec_data              = mean_subset_data(absSWRinRec_idx,:);
% awake_swr_rate_per_min_data   = mean_subset_data(awake_swr_rate_per_min_idx,:);
% propClustSWR_data             = mean_subset_data(propClustSWR_idx,:);
% swr_dur_in_sec_data           = mean_subset_data(swr_dur_in_sec_idx,:);
% swrAmp_data                   = mean_subset_data(swrAmp_idx,:);
% nr_swr_clusters_in_rec_data   = mean_subset_data(nr_swr_clusters_in_rec,:);
% prop_swr_clusters_in_rec_data = mean_subset_data(prop_swr_clusters_in_rec,:);
% ISI_data                      = mean_subset_data(ISI,:);
% 
% % Convert from cell to array
% absSWRinRec_array              = cell2mat(absSWRinRec_data);
% awake_swr_rate_per_min_array   = cell2mat(awake_swr_rate_per_min_data);
% propClustSWR_array             = cell2mat(propClustSWR_data);
% swr_dur_in_sec_array           = cell2mat(swr_dur_in_sec_data);
% swrAmp_array                   = cell2mat(swrAmp_data);
% nr_swr_clusters_in_rec_array   = cell2mat(nr_swr_clusters_in_rec_data);
% prop_swr_clusters_in_rec_array = cell2mat(prop_swr_clusters_in_rec_data);
% ISI_array                      = cell2mat(ISI_data);
% 
% mean_feat_awake = {absSWRinRec_array, awake_swr_rate_per_min_array, propClustSWR_array,...
%             swr_dur_in_sec_array, swrAmp_array, nr_swr_clusters_in_rec_array,...
%             prop_swr_clusters_in_rec_array,ISI_array};
