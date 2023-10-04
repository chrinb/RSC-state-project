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

imaging_sampling_rate = find_imaging_framerate(sData);
n_rois                = size(dff, 1);
sig_transients        = sData.imdata.roiSignals(2).([params.cell_type, '_sig_transients']);


% Preallocate
[transient_duration_sec, transient_times, transient_amplitudes] = deal( cell(size(dff,1),1));

% Loop over ROIs
for roi_nr = 1:n_rois
    
    transient     = sig_transients(roi_nr,:);
    transient_log = sData.analysis.transients.sig_transients_logmat(roi_nr,:);

    % Compute transient durations
    [eventStartIdx, eventStopIdx ] = findTransitions( transient_log);
    transient_duration_sec{roi_nr} = (eventStopIdx - eventStartIdx)./imaging_sampling_rate;

    % Find peak DF/F amplitude per transient
    transient_amplitudes{roi_nr, 1} = cell2mat( arrayfun(@(x) max(transient(eventStartIdx(x):eventStopIdx(x))), (1:size(eventStartIdx,2)), 'uni', 0));

    % Find time points of transient onsets
    transient_times{roi_nr,1} = eventStartIdx;
end

% Store in sData
sData.analysis.transients.([params.cell_type, '_transient_dur_sec']) = transient_duration_sec;
sData.analysis.transients.([params.cell_type, '_transient_ampl'])    = transient_amplitudes;
sData.analysis.transients.([params.cell_type, '_transient_times'])   = transient_times;
