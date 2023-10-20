function output = find_across_day_ROI_idx(across_session_idx, roi_idx)

% Written by Christoffer Berge | Vervaeke lab

%{
When imaging the same FOV repeatedly across sessions, some cells will appear 
in multiple sessions while others only in one. This function 
(1) determines which ROIs are present on multiple sessions
(2) finds which of those sessions the ROI is present in
(3) finds the index of that ROI in each session
%}

across_session_idx   = cell2mat( across_session_idx);
unique_FOV_list      = unique( across_session_idx);
total_nr_of_sessions = 1:size(across_session_idx,1);

% Loop over all unique FOVs
for unique_FOV_nr = 1:numel(unique_FOV_list)

    sessions_to_merge_idx = across_session_idx == unique_FOV_nr;

    if sum(sessions_to_merge_idx)  > 1
        roi_indices           = roi_idx(sessions_to_merge_idx);
    
        % Indices of ROIs to analyze
        temp = horzcat(roi_indices{:});
    
        % Find duplicate indices (to be be averaged)
        min_val = min(temp);
        max_val = max(temp);
    
        idx_of_duplicates = zeros(1,  size(min_val:max_val,2));
        for j = min_val:max_val
            if sum(temp == j) > 1
                idx_of_duplicates(1, j) = j;
            end
            idx_of_duplicates(idx_of_duplicates == 0) = [];
        end
    
        % Find location of these indices
        session_to_merge = total_nr_of_sessions(sessions_to_merge_idx);
        counter = 1;
        for k = idx_of_duplicates
       
            for session_nr = 1:size(session_to_merge,2)
                if isempty(find(roi_indices{session_nr, 1}==k,1))
                    rois_to_average(counter, session_nr) = NaN;
                else
                    rois_to_average(counter, session_nr) = find(roi_indices{session_nr, 1}==k);
                end
            end
            rois_to_average(counter,session_nr+1) = k;
            counter = counter + 1;
        end
    
        for session_idx = 1:size(session_to_merge,2)
            session_nr = session_to_merge(session_idx);
            all_table_text{session_idx} = ['ROI index in session ', num2str(session_nr)];
        end
    
        all_table_text{1, session_idx+1} = 'ROI nr';
        table_of_indices = array2table(rois_to_average, 'VariableNames',all_table_text);
    
        output{unique_FOV_nr} = table_of_indices;
        clear rois_to_average
    end
end


