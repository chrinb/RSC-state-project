function params = select_roi_im_epochs(sData, params)

%{
Select epochs for creating mean ROI images in a recording. 
%}

% Select epochs either as before-during-after REM episode, or as
% beginning and end of recording
if strcmp(params.rem, 'yes')
   
    % Automatic selection of before-during-after REM epoch if not specified
    % by user           
    frames  = sData.daqdata.frame_onset_reference_frame;
    if sum( isnan(params.EpochFrameStart)) == 3
        rem_start_stop = frames( rem_sleep(sData));

        % if more than 1 REM ep, choose whichever clostes to the middle of recording as the
        % during-epoch and avoid the other(s) in the before-after epochs
        if numel(rem_start_stop(:,1)) > 1
            n_frames = sData.daqdata.frame_onset_reference_frame(end);
            [~, idx]    = min(abs(rem_start_stop(:,1)-n_frames));
            
            if idx == 2
                before_start_idx = rem_start_stop(1,2);
                after_end_idx = frames(end);
            elseif idx == 1
                before_start_idx = 1;
                after_end_idx = rem_start_stop(2,1);
            end
            % [~, max_id] = max( rem_start_stop(:,2)-rem_start_stop(:,1) );
        else
            idx = 1;
            before_start_idx = 1;
            after_end_idx = frames(end);
        end

        rem_during     = median( rem_start_stop(idx,:) );
        rem_before     = median( [before_start_idx, rem_start_stop(idx,1)] );
        rem_after      = median( [rem_start_stop(idx, end), after_end_idx]);
    
        params.EpochFrameStart = round([rem_before, rem_during, rem_after]);
    end

elseif strcmp(params.rem, 'no')

    params.EpochFrameStart = [1 17000];
end