function  [data, data_zscore, data_SE, data_zscore_SE] = get_mean_event_data(varargin)

% Written by Christoffer Berge

% Calculate mean activity during SWRs and spindles. First input argument
% has to be a struct containing the mean data (e.g., dF/F or deconvolved
% dF/F) from several sessions. Second input argument is the name of the
% struct field this function will search for. For example, to find all mean
% SWR spindle-coupled sessions type "coupled". If sessions from multiple
% mice are loaded in struct, first specify the mouse name ("m6120"), then
% specify SWR type as third input ("awake"). 

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
data_idx        = (1:2:length(mean_subset_data));
data_zscore_idx = (2:2:length(mean_subset_data));

data        = mean_subset_data(data_idx,:);
data_zscore = mean_subset_data(data_zscore_idx,:);

% Convert from cell to array
data        = cell2mat(data);
data_zscore = cell2mat(data_zscore);

% Calculate standard error of the mean
data_SE        = std(data, 'omitnan') ./ sqrt(size(data,1));
data_zscore_SE = std(data_zscore, 'omitnan') ./ sqrt(size(data_zscore,1));
