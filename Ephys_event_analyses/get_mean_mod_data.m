function output = get_mean_mod_data(varargin)

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
ROIs_activated_idx         = (1:8:length(mean_subset_data));
ROIs_activated_idx_idx     = (2:8:length(mean_subset_data));
ROIs_suppressed_idx        = (3:8:length(mean_subset_data));
ROIs_suppressed_idx_idx    = (4:8:length(mean_subset_data));
activated_sig_to_plot_idx  = (5:8:length(mean_subset_data));
suppressed_sig_to_plot_idx = (6:8:length(mean_subset_data));
ROIs_unclassified_idx      = (7:8:length(mean_subset_data));
ROIs_unclassified_idx_idx  = (8:8:length(mean_subset_data));


ROIs_activated_idx_data         = mean_subset_data(ROIs_activated_idx,:);
ROIs_activated_idx_idx_data     = mean_subset_data(ROIs_activated_idx_idx,:);
ROIs_suppressed_idx_data        = mean_subset_data(ROIs_suppressed_idx,:);
ROIs_suppressed_idx_idx_data    = mean_subset_data(ROIs_suppressed_idx_idx,:);
activated_sig_to_plot_idx_data  = mean_subset_data(activated_sig_to_plot_idx,:);
suppressed_sig_to_plot_idx_data = mean_subset_data(suppressed_sig_to_plot_idx,:);
ROIs_unclassified_idx_data      = mean_subset_data(ROIs_unclassified_idx,:);
ROIs_unclassified_idx_idx_data  = mean_subset_data(ROIs_unclassified_idx_idx,:);

% Convert from cell to array
ROIs_activated_idx_array         = cell2mat(ROIs_activated_idx_data);
ROIs_activated_idx_idx_array     = cell2mat(ROIs_activated_idx_idx_data);
ROIs_suppressed_idx_array        = cell2mat(ROIs_suppressed_idx_data);
ROIs_suppressed_idx_idx_array    = cell2mat(ROIs_suppressed_idx_idx_data);
activated_sig_to_plot_idx_array  = cell2mat(activated_sig_to_plot_idx_data);
suppressed_sig_to_plot_idx_array = cell2mat(suppressed_sig_to_plot_idx_data);
ROIs_unclassified_idx_array      = cell2mat(ROIs_unclassified_idx_data);
ROIs_unclassified_idx_idx_array  = cell2mat(ROIs_unclassified_idx_idx_data);


output = {ROIs_activated_idx_array, ROIs_activated_idx_idx_array, ROIs_suppressed_idx_array,...
            ROIs_suppressed_idx_idx_array, activated_sig_to_plot_idx_array,suppressed_sig_to_plot_idx_array,...
            ROIs_unclassified_idx_array,ROIs_unclassified_idx_idx_array};
