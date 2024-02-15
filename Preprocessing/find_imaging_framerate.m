function imaging_sampling_rate = find_imaging_framerate(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Function that finds the two photon imaging frame rate

sData = varargin{1,1};

if isfield(sData, 'daqdataB') && isfield(sData, 'behavior')
    imaging_sampling_rate = sData.behavior.meta.imagingSamplingRate;
% elseif isfield(sData.behavior.meta, 'imagingSamplingRate')
%     imaging_sampling_rate = sData.behavior.meta.imagingSamplingRate;

elseif isfield(sData, 'daqdata')
    frame_signal = sData.daqdata.frameSignal;

    

    if ~isempty(frame_signal)
        frame_start           = diff(frame_signal)== 32;     
        frame_start_index     = find(frame_start);
    else
        frame_start_index = find(sData.daqdata.frame_onset);
    end

    imaging_sampling_rate = 1 / ( (mean(diff(frame_start_index))) / sData.daqdata.metadataFromLV.Sampling_rate_downsampled  ); 
end
% If piezo recording divide by 4
if nargin > 1
    imaging_sampling_rate = imaging_sampling_rate/4;
end
