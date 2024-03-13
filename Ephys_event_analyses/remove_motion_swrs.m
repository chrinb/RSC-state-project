function swr_keep_idx = remove_motion_swrs(sData)

%{
Remove SWRs occurring during locomotion bouts as identified via paw
movements for the PRE- and POST-RUN sessions.
%}

% Remove SWRs occurring after 2P imaging stop
dff                   = sData.imdata.roiSignals(2).newdff;
swr_in_2p_idx         = sData.ephysdata.frameRipIdx < size(dff,2);
frame_swr_idx_trimmed = sData.ephysdata.frameRipIdx(swr_in_2p_idx);

% Check whether SWR timestaps (in 2P frame time) occur within the QW
% epochs of the movement threshold vector and keep SWRs that do
rec_length_vector = 1:length(sData.analysis.movement_threshold);
samples_to_keep   = rec_length_vector(~sData.analysis.movement_threshold);
swr_keep_idx      = ismember(frame_swr_idx_trimmed, samples_to_keep);
