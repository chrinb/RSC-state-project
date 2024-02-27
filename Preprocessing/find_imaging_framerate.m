function srate = find_imaging_framerate(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Function that finds the two photon imaging frame rate

sData = varargin{1,1};

if isfield(sData, 'daqdataB') && isfield(sData, 'behavior')
    srate = sData.behavior.meta.imagingSamplingRate;
% elseif isfield(sData.behavior.meta, 'imagingSamplingRate')
%     imaging_sampling_rate = sData.behavior.meta.imagingSamplingRate;

elseif isfield(sData, 'daqdata')
    frame_signal = sData.daqdata.frameSignal;
    

    if ~isempty(frame_signal)
        frame_start           = diff(frame_signal)== 32;     
        frame_start_index     = find(frame_start);
        
        % check for missing frames
        [vals, val_idx] = unique(diff(frame_start_index));

        med_val = (median(diff(frame_start_index)));

        missing_frame_idx = [];
        for i = 1:numel(vals)
            if vals(i) > med_val + 2
                missing_frame_idx(i) = i;
            end
        end
        missing_frame_idx(missing_frame_idx==0) = [];

        if numel(missing_frame_idx) > 1
            errordlg( ['Several indices of missing frames for ', sData.sessionInfo.sessionID,', investigate'] )
            
        elseif numel(missing_frame_idx)== 1
            msgbox( ['Missing frames at index ', num2str(val_idx(missing_frame_idx)), ' for ', sData.sessionInfo.sessionID,' will remove subsequent frames'] )
            frame_start_index(val_idx(missing_frame_idx):end) = [];
        end

    else
        frame_start_index = find(sData.daqdata.frame_onset);
    end

    srate = 1 / ( (mean(diff(frame_start_index))) / sData.daqdata.metadataFromLV.Sampling_rate_downsampled  ); 
    msgbox( ['Sample rate for ', sData.sessionInfo.sessionID,' is ', num2str(srate)])
end
% If piezo recording divide by 4
if nargin > 1
    srate = srate/4;
end
