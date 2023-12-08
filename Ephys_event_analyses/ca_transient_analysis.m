function sData = ca_transient_analysis(sData, params)

% Written by Christoffer Berge | Vervaeke lab

%{
Function that calculates duration, amplitude, and onset of significant
calcium transients
%}

% Also include "integrated DF/F per transient"??

[pc_rois, in_rois] = remove_cells_longitudinal(sData, params);

% Select data
switch params.cell_type
    case 'axon'
    sig_transients = sData.imdata.roiSignals(2).mergedAxons_sig_transients;
    txt = 'axon';
    case 'pc'
    sig_transients = sData.imdata.roiSignals(2).pc_sig_transients;
    txt = 'pc';
    case 'in'
    sig_transients = sData.imdata.roiSignals(2).in_sig_transients;
    txt = 'in';
    % case 'pc'
    % sig_transients = sData.imdata.roiSignals(2).all_sig_transients(pc_rois,:);
    % txt = 'pc';
    % case 'in'
    % sig_transients = sData.imdata.roiSignals(2).all_sig_transients(in_rois,:);
    % txt = 'in';
%     casre 'all'
%     dff = sData.imdata.roiSignals(2).newdff;
end

imaging_sampling_rate = find_imaging_framerate(sData);
n_rois                = size(sig_transients, 1);
% sig_transients        = sData.imdata.roiSignals(2).([params.cell_type, '_sig_transients']);


% Preallocate
[transient_duration_sec, transient_times, transient_amplitudes] = deal( cell(size(sig_transients,1),1));

% Loop over ROIs
for roi_nr = 1:n_rois
    
    transient     = sig_transients(roi_nr,:);
    transient_log = sData.analysis.transients.([txt, '_sig_transients_logmat'])(roi_nr,:);

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
