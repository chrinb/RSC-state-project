function sData = ca_transient_analysis(sData, params)

% Written by Christoffer Berge | Vervaeke lab

% Function that for each ROI DF/F trace with significant Ca2+ transients
% computes: (1) transient peak DF/F, (2) transients per second, (3)
% transient duration (in seconds), (4) time points of transient onsets.

% Also include "integrated DF/F per transient"??

% TO DO: is it worth splitting excitatory and inhibitory ROIS??? will need
% to save them separately...


%% Get exctitatory and inhibitory indices
[pc_rois, in_rois] = remove_cells_longitudinal(sData);

% Select data
switch params.cell_type
    case 'axon'
     dff    = sData.imdata.roiSignals(2).mergedAxonsDffFilt;
    case 'in'
    dff = sData.imdata.roiSignals(2).newdff(in_rois,:);
    case 'pc'
    dff = sData.imdata.roiSignals(2).newdff(pc_rois,:);
    case 'all'
    dff = sData.imdata.roiSignals(2).newdff;
end

n_rois         = size(dff,1);
sig_transients = sData.analysis.transients.roi_signal_sig_transients;
rec_length     = size(dff, 2)/31;

% Preallocate
durations_in_sec    = cell(size(dff,1),1);
% all_durations       = [];
all_amplitudes      = cell(size(dff,1),1);
all_peak_amplitudes = [];
transients_per_sec  = zeros(size(dff,1),1);
transient_times     = cell(size(dff,1),1);

% Loop over ROIs
for roi_nr = 1:n_rois
    
    transient = sig_transients(roi_nr,:);

    % Compute transient durations
    temp                            = ~isnan( transient);
    [eventStartIdx, eventStopIdx ]  = findTransitions( double( temp));
    durations_in_sec{roi_nr,1}      = ( (eventStopIdx+1) - eventStartIdx)./31;
%     all_durations                   = [all_durations, durations_in_sec{roi_nr,1} ];

    % Find peak DF/F amplitudes
    
    % Loop over transients
    temp_ampl = [];
    for transient_nr = 1:size(eventStartIdx,2)
        
        % Find max value of current significant transient
        max_dff_per_transient     = max( transient( eventStartIdx( transient_nr):eventStopIdx( transient_nr)));
%         all_amplitudes{roi_nr, transient_nr} = max_dff_per_transient;
        temp_ampl = [temp_ampl; max_dff_per_transient ];
%         all_peak_amplitudes       = [all_peak_amplitudes, max_dff_per_transient];
    end
    all_amplitudes{roi_nr, 1} = temp_ampl;


    % Compute transients per second
    
%     total_transient_n          = size( eventStartIdx,2);
%     transients_per_sec(roi_nr) = total_transient_n/rec_length;

    % Find time points of transient onsets
    transient_times{roi_nr,1} = eventStartIdx;
end

sData.analysis.transients.durations_in_sec    = durations_in_sec;
% sData.analysis.transients.all_durations       = all_durations;
sData.analysis.transients.all_amplitudes      = all_amplitudes;
% sData.analysis.transients.all_peak_amplitudes = all_peak_amplitudes;
% sData.analysis.transients.transients_per_sec  = transients_per_sec;
sData.analysis.transients.transient_times     = transient_times;
