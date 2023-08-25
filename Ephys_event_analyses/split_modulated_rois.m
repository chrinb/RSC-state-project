function sorted_rois = split_modulated_rois(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that takes event-modulated ROIs from "event_mod_an" and (1) split
% them into activated or suppressed, and (2) get the corresponding DF/F
% signals from the same ROIs

checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

if size(varargin) == [1,1]
    varargin = varargin{1,1};
end

All_ROIs                          = varargin{1,1};
win_length                        = varargin{1,2};
signal_event_activity_mean        = varargin{1,3};
signal_event_activity_mean_zscore = varargin{1,4};
opts                              = varargin{1,5};
sData                             = varargin{1,6};
event_idx                         = varargin{1,7};
nr_of_seconds                     = varargin{1,8};


%% Find corresponding DF/F signals 

% Depending on which signal was used for modulation analysis, extract the
% other signal (e.g., if DF/F extract deconv. and vice versa)
if checkParameter(opts.signal_type, 1 , 'deconv')
    opts.signal_type = 'dff';
elseif checkParameter(opts.signal_type, 2 , 'dff')
    opts.signal_type = 'deconv';
end


% Get ROI signals
[signal, text, opts, label3, rois_for_an, ~, cells_to_exclude] = ...
get_roi_signals_from_sData(sData, opts );


% Preallocate again, this time for DF/F signals
signal_event_activity                 = zeros(length(event_idx), win_length);
signal_event_activity_zscore          = zeros(length(event_idx), win_length);
signal_event_activity_mean_dff        = zeros(size(rois_for_an,2), win_length);
signal_event_activity_mean_dff_zscore = zeros(size(rois_for_an,2), win_length);

for roinr = rois_for_an
    roi_signal        = signal(roinr, :);
    roi_signal_zscore = zscore(roi_signal);
    
    [signal_event_activity,signal_event_activity_zscore] = ...
        extract_avg_activity(event_idx, roi_signal, roi_signal_zscore, signal_event_activity,...
        signal_event_activity_zscore,nr_of_seconds, opts);
    
    % Find mean/median SWR-aligned activity for each ROI
    signal_event_activity_mean_dff(roinr,:)        = mean(signal_event_activity, 'omitnan');
    signal_event_activity_mean_dff_zscore(roinr,:) = mean(signal_event_activity_zscore, 'omitnan');
       
end

%% Find activated and suppressed ROIs

all_modulated_ROIs = nonzeros(All_ROIs);

% loop over modulated ROIs and classify as activated or suppressed based on
% comparing mean activity in a baseline window (-3 to -2 sec) vs a test
% window (-500ms to +500ms)
base_win       = (1:win_length);
frames_in_1sec = 31;
test_window    = select_mod_win(win_length, opts.event_type);
% test_window    = (win_length/2) - frames_in_1sec/2 : (win_length/2) + frames_in_1sec/2;

[ROIs_activated_dff, ROIs_activated_dff_zscore, ROIs_activated_idx,...
    ROIs_suppressed_dff, ROIs_suppressed_dff_zscore, ROIs_suppressed_idx,...
    ROIs_unclassified, ROIs_unclassified_zscore, ROIs_unclassified_idx] = deal([]);

% Loop over modulated ROIs
for neuron_n = all_modulated_ROIs'
    
    % Check if mean activity in baseline window is greater or smaller than
    % mean activity in user defined window, and assign ROIs as suppressed
    % or activated, respectively. Remaining modulated ROIs are
    % unclassified.

    % TO DO difference have to be of some meaningful magnitude. slight
    % difference will lead to weird results.
    if mean( signal_event_activity_mean_dff(neuron_n,base_win)) > mean( signal_event_activity_mean_dff(neuron_n,test_window))
        ROIs_s                     = signal_event_activity_mean_dff( neuron_n,:);
        ROIs_s_zscore              = signal_event_activity_mean_dff_zscore( neuron_n,:);
        ROIs_suppressed_dff        = [ROIs_suppressed_dff; ROIs_s];
        ROIs_suppressed_dff_zscore = [ROIs_suppressed_dff_zscore; ROIs_s_zscore];
        ROIs_suppressed_idx        = [ROIs_suppressed_idx; neuron_n];
    
    elseif mean( signal_event_activity_mean_dff(neuron_n,base_win)) < mean( signal_event_activity_mean_dff(neuron_n, test_window))
        ROIs_a                    = signal_event_activity_mean_dff( neuron_n,:);
        ROIs_a_zscore             = signal_event_activity_mean_dff_zscore( neuron_n,:);
        ROIs_activated_dff        = [ROIs_activated_dff; ROIs_a];
        ROIs_activated_dff_zscore = [ROIs_activated_dff_zscore; ROIs_a_zscore];
        ROIs_activated_idx        = [ROIs_activated_idx; neuron_n];
    
%     else
%         ROIs_u                   = signal_swr_activity_mean( neuron_n,:);
%         ROIs_u_zscore            = signal_swr_activity_mean_zscore( neuron_n, :);
%         ROIs_unclassified        = [ROIs_unclassified; ROIs_u];
%         ROIs_unclassified_zscore = [ROIs_unclassified_zscore; ROIs_u_zscore];
%         ROIs_unclassified_idx    = [ROIs_unclassified_idx; neuron_n];
    end
end
    

% Get deconvolved signal from modulated ROIs
ROIs_suppressed_dec          = signal_event_activity_mean(ROIs_suppressed_idx,:);
ROIs_suppressed_dec_zscore   = signal_event_activity_mean_zscore(ROIs_suppressed_idx,:);
ROIs_activated_dec           = signal_event_activity_mean(ROIs_activated_idx,:);
ROIs_activated_dec_zscore    = signal_event_activity_mean_zscore(ROIs_activated_idx,:);
% ROIs_unclassified_dff        = signal_swr_activity_mean(ROIs_activated_idx,:);
% ROIs_unclassified_dff_zscore = signal_swr_activity_mean_zscore(ROIs_activated_idx,:);

sorted_rois = {ROIs_activated_dff,    ROIs_activated_dff_zscore,    ROIs_activated_dec,    ROIs_activated_dec_zscore,    ROIs_activated_idx...
               ROIs_suppressed_dff,   ROIs_suppressed_dff_zscore,   ROIs_suppressed_dec,   ROIs_suppressed_dec_zscore,   ROIs_suppressed_idx};
%                ROIs_unclassified, ROIs_unclassified_zscore, ROIs_unclassified_dff, ROIs_unclassified_dff_zscore, ROIs_unclassified_idx};