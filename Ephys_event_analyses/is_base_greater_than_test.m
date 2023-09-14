function check = is_base_greater_than_test(signal_event_activity_mean_dff, roi_nr, base_win, test_window)

% Written by Christoffer Berge | Vervaeke lab

% Check if mean activity in baseline window is greater than test window

check = mean( signal_event_activity_mean_dff(roi_nr, base_win)) > mean( signal_event_activity_mean_dff(roi_nr, test_window));