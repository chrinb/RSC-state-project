function output = concatenate_sleep_data(varargin)

n_sessions = size(varargin,2);
n_rois     = size(varargin{1, 1}.sort_idx,1);
% Loop over ROIs
for roi_nr = 1:n_rois
    
    [data_unclassified, data_single, data_coupled] = deal([]);
    [swr1, swr2, swr3] = deal(0);
    % Loop over sessions
    for session_nr = 1:n_sessions
    
        % Get ROI data from current session
        temp_unclassified = varargin{1, session_nr}.all_data_zscored_unclassified(roi_nr,:,:);
        temp_single       = varargin{1, session_nr}.all_data_zscored_single(roi_nr,:,:);
        temp_coupled      = varargin{1, session_nr}.all_data_zscored_coupled(roi_nr,:,:);
        
        % Count nr of SWRs
        temp_swr1 = size(temp_unclassified,2);
        temp_swr2 = size(temp_single,2);
        temp_swr3 = size(temp_coupled,2);

        % Concatenate ROI data across sessions
        data_unclassified = [data_unclassified, temp_unclassified];
        data_single       = [data_single, temp_single];
        data_coupled      = [data_coupled, temp_coupled];

        swr1 = swr1 + temp_swr1;
        swr2 = swr2 + temp_swr2;
        swr3 = swr3 + temp_swr3;
    end
    
    concatenated_unclassified_data(roi_nr,:,:) = mean(data_unclassified, 2);
    concatenated_single_data(roi_nr,:,:)       = mean(data_single, 2);
    concatenated_coupled_data(roi_nr,:,:)      = mean(data_coupled, 2);
end


output.signal_swr_activity_mean_unclassified_zscore = concatenated_unclassified_data;
output.signal_swr_activity_mean_single_zscore       = concatenated_single_data;
output.signal_swr_activity_mean_coupled_zscore      = concatenated_coupled_data;
output.swr_nr       = [swr1, swr2, swr3]; 
