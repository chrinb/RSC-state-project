function [ROIs_activated_idx, ROIs_suppressed_idx] = split_modulated_rois(signal_event_activity_mean_zscore, win_length, params, event_mod_ROI_idx)

% Written by Christoffer Berge | Vervaeke Lab

% Find positively or negatively ROIs during some event by comparing mean
% peri-event activity to mean baseline activity

% Define baseline and peri-event (test) wndow
base_win       = (1:31*2);
test_window    = select_mod_win(win_length, params.event_type);

[ROIs_activated_idx, ROIs_suppressed_idx] = deal([]);

% Loop over modulated ROIs
for roi_nr = event_mod_ROI_idx

    if is_base_greater_than_test(signal_event_activity_mean_zscore, roi_nr, base_win, test_window)
        ROIs_suppressed_idx        = [ROIs_suppressed_idx; roi_nr];
    else
        ROIs_activated_idx        = [ROIs_activated_idx; roi_nr];
    end
end

