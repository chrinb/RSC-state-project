function peri_event_activity = swr_aligned_ephys(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Calculate peri-SWR activity for various ephys signals. 

sData     = varargin{1,1};
params    = varargin{1,2};
sleep_idx = varargin{1,3};

% Load and assign variables
nr_of_seconds  = 1;
win_length     = (nr_of_seconds*2500*2)+1;
select_swr_idx = sData.ephysdata.absRipIdx; 

% Select signal
if strcmp(params.signal_type, 'ECoG')
    signal = sData.ephysdata2.lfp';
elseif strcmp(params.signal_type, 'EMG')
    signal = sData.ephysdata3.lfp';
elseif strcmp(params.signal_type,'Running speed')
    signal = sData.daqdata.runSpeed';
elseif strcmp(params.signal_type,'LFP')
    signal = sData.ephysdata.lfp';
end

% Z-score signal (not running speed)
if strcmp(params.zscore, 'yes') && ~strcmp(params.signal_type, 'Running speed')
    % signal = zscore(signal);
    % signal = signal-mean(signal);

end

%% Select SWRs for analysis

% If sleep session, merge all NREM SWRs (spindle-coupled and
% spindle-uncoupled)
if sleep_idx == 1
    select_swr_idx = sort([sData.ephysdata.NREM_spindle_uncoupled_swr, sData.ephysdata.spindle_coupled_swr]);
end

% Get the indicies of user specified SWR types
event_idx = get_swr_idx(params.swr_for_analysis, sData, select_swr_idx, params);

n_ripples = length(event_idx);

%% Run main analysis

% Preallocate
peri_event_activity = deal( zeros(n_ripples, win_length));

params.baseSub = 'subtract';
% Loop over nr of SWRs
for swr_nr = 1:n_ripples
    
    % Extract peri-SWR ephys data
    peri_event_activity = extract_avg_activity(event_idx, signal, peri_event_activity, nr_of_seconds, params);
end
 
