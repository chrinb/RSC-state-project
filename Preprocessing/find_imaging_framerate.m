function imaging_sampling_rate = find_imaging_framerate(sData)

% Written by Christoffer Berge || Vervaeke lab

% Function that finds the two photon imaging frame rate

if ~isempty(sData.daqdata.frameSignal)
    frame_start           = diff(sData.daqdata.frameSignal)== 32;     
    frame_start_index     = find(frame_start);
else
    frame_start_index = find(sData.daqdata.frame_onset);
end
imaging_sampling_rate = 1 / ( (mean(diff(frame_start_index))) / sData.daqdata.metadataFromLV.Sampling_rate_downsampled  ); 