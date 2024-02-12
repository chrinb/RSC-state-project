function peri_event_activity = extract_avg_activity(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that extracts average ROI activity during a particular event
% (SWRs, sleep spindles, slow waves)


% Get variables 
event_idx                    = varargin{1,1}; % SWR idx, spindle idx, slow wave idx, etc.
signal                       = varargin{1,2}; 
% roi_signal_zscore              = varargin{1,3};
peri_event_activity        = varargin{1,3};
% signal_event_activity_zscore = varargin{1,5};
nr_of_seconds                = varargin{1,4};
params                       = varargin{1,5};

% anon function for parameter check
checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Define baseline window. For two-photon fluorescence data window
% is defines as the mean of the 31 (1 sec) first frames of the
% peri-event window. For ephys, its the first 2500 samples (1 sec )
% of the peri-event window

if strcmp(params.signal_type, 'EMG') || strcmp(params.signal_type, 'ECoG') || strcmp(params.signal_type, 'Running speed')
     baseline_win = 1:2500;
     one_sec = 2500;
else
    baseline_win = 1:31;
    one_sec = 31;
end

% Loop over nr of events
for event_nr = 1:length(event_idx)

    % find event window start end
    event_window_start = event_idx(event_nr) - (nr_of_seconds*one_sec); 
    event_window_end   = event_idx(event_nr) + (nr_of_seconds*one_sec); 

    % Extract snippets of data around event. Skip events at the beginning 
    % or end with a window shorter than length specified in seconds by user above
    if event_window_start > 1 && event_window_end < size(signal,2)
       
        % Data
        peri_event_activity(event_nr, :) = ...
            signal(event_window_start:event_window_end);
        % Z-scored data
%         signal_event_activity_zscore(event_nr, :) = ...
%             roi_signal_zscore(event_window_start:event_window_end); 
    end
    
    % Baseline subtraction
    if checkParameter(params.baseSub, 1, 'subtract')
        baselineDFF                              = mean(peri_event_activity(event_nr, baseline_win),'omitnan');
%         baselineDFF_zscore                       = mean(signal_event_activity_zscore(event_nr, baseline_win), 'omitnan');
        peri_event_activity(event_nr,:)        = peri_event_activity(event_nr,:)-baselineDFF;
%         signal_event_activity_zscore(event_nr,:) = signal_event_activity_zscore(event_nr,:)-baselineDFF_zscore;
    end
end
