function mean_feat_sleep = get_mean_swr_feat_sleep(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that extracts data from a multi-session struct containing different
% variables that repeats across sessions and pools together the repeating 
% variable types. For example, it locates every struct field in the multi-session
% struct that stores the nr of SWRs in each session data and then groups
% these values into  

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

% Initialize data matrix
% mean_subset_data = [];
% 
% % Loop over field names 
% for i=1:length(f)
%     % if two field names are specified, check that each field name contains
%     % both strings
%     if size(varargin,2) == 3
%         split_field = split(f{i}, "_");
%         if contains( split_field(1), fieldName) && contains( split_field(5), fieldName2)
%             mean_subset_data = vertcat(mean_subset_data, data_cell(i) );
%         end
% 
%     elseif contains( f{i}, strtrim(fieldName),'IgnoreCase',true)
%         mean_subset_data = vertcat(mean_subset_data, data_cell(i) );
%     end
% end

% Create a new cell array where first row contains unique variable names
mean_feat_sleep = cellstr(unique_variable_names');

% Loop over nr of unique variable names 
for i = 1:nr_of_unique_field_names

    % The data in the multi-session struct repeats every 
    % "nr-of-unique-variable-names"-th time. Use this nr to index repeating
    % variables. 
    data_idx             = i:nr_of_unique_field_names:length(mean_subset_data);
    data_cell            = mean_subset_data(data_idx,:);
    data_arr             = cell2mat(data_cell);
    mean_feat_sleep{2,i} = data_arr;
end

% output = mean_subset_data;
% Index data
% absSWRinRec_idx                  = (1:15:length(mean_subset_data));
% propNREMClustSWR_idx             = (2:15:length(mean_subset_data));
% propClustSWR_idx                 = (3:15:length(mean_subset_data));
% NREMswrDur_idx                   = (4:15:length(mean_subset_data));
% NREMswrRate_idx                  = (5:15:length(mean_subset_data));
% NREMswrAmp_idx                   = (6:15:length(mean_subset_data));
% nrSpindleCoupledSWR_idx          = (7:15:length(mean_subset_data));
% propSpindleCoupledSWR_idx        = (8:15:length(mean_subset_data));
% nr_swr_clusters_in_rec           = (9:15:length(mean_subset_data));
% prop_swr_clusters_in_rec         = (10:15:length(mean_subset_data));
% nr_swr_clusters_in_nrem          = (11:15:length(mean_subset_data));
% prop_swr_cluster_in_nrem         = (12:15:length(mean_subset_data));
% nr_cluster_swr_in_sleep_spindles = (13:15:length(mean_subset_data));
% prop_swr_cluster_in_spindle      = (14:15:length(mean_subset_data));
% ISI                              = (15:15:length(mean_subset_data));
% 
% absSWRinRec_data      = mean_subset_data(absSWRinRec_idx,:);
% propNREMClustSWR_data = mean_subset_data(propNREMClustSWR_idx,:);
% propClustSWR_data     = mean_subset_data(propClustSWR_idx,:);
% 
% NREMswrDur_data  = mean_subset_data(NREMswrDur_idx,:);
% NREMswrRate_data = mean_subset_data(NREMswrRate_idx,:);
% NREMswrAmp_data  = mean_subset_data(NREMswrAmp_idx,:);
% 
% nrSpindleCoupledSWR_data   = mean_subset_data(nrSpindleCoupledSWR_idx,:);
% propSpindleCoupledSWR_data = mean_subset_data(propSpindleCoupledSWR_idx,:);
% % ISI_data = mean_subset_data(ISI,:)
% % Convert from cell to array
% absSWRinRec_array      = cell2mat(absSWRinRec_data);
% propNREMClustSWR_array = cell2mat(propNREMClustSWR_data);
% propClustSWR_array     = cell2mat(propClustSWR_data);
% 
% NREMswrDur_array  = cell2mat(NREMswrDur_data);
% NREMswrRate_array = cell2mat(NREMswrRate_data);
% NREMswrAmp_array  = cell2mat(NREMswrAmp_data);
% 
% nrSpindleCoupledSWR_array   = cell2mat(nrSpindleCoupledSWR_data);
% propSpindleCoupledSWR_array = cell2mat(propSpindleCoupledSWR_data);
% 
% 
% output = {absSWRinRec_array, propNREMClustSWR_array, propClustSWR_array,...
%             NREMswrDur_array, NREMswrRate_array,NREMswrAmp_array,...
%             nrSpindleCoupledSWR_array,propSpindleCoupledSWR_array};
