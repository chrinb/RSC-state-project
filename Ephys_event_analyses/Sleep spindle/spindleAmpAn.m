function spindleAmpAn(varargin)

sData = varargin{1,1};

% if isempty(swr_idx)
%     swr_idx = true([1,length(sData.ephysdata.absRipIdx)]);
% else
%     swr_idx = swr_idx; % for specifying subset of SWRs, for example NREM SWRs.
% end
i = 1;
ripWindow = 1; %define how many seconds before and after ripple you want to look at 

% enumerate_swr = 1:length(sData.ephysdata.absRipIdx);
% swr_idx_to_use = enumerate_swr(swr_idx);
% for i = 1:enumerate_swr
    
spindles = sData.ephysdata2.NREMspindleSnips;
% spinLocs = sData.ephysdata2.frameRipIdx;
imgFs = 31;
nFrames = max(sData.daqdata.frame_onset_reference_frame);
nFrames_window = imgFs*ripWindow;

fs = 2500;
nyquistFs = fs/2;
filter_kernel = fir1(600,[8 18]./nyquistFs);
time = linspace(0,1,2500);


dF = sData.imdata.roiSignals(2).newdff;

%eliminate all ripples that don't have enough time before or after to
%perform analysis (happen too early or too late in the recording)
minTimeWindow = nFrames_window + 1;
maxTimeWindow = nFrames-(nFrames_window);
% outOfBoundsIdx = find(ripLocs < minTimeWindow);
% ripLocs(outOfBoundsIdx) = [];
% ripples(outOfBoundsIdx) = [];
% outOfBoundsIdx = find(ripLocs > maxTimeWindow);
% ripLocs(outOfBoundsIdx) = [];
% ripples(outOfBoundsIdx) = [];

for j = 1:length(spindles)
        lfp = spindles(j).lfp;
        filtered_lfp = filtfilt(filter_kernel,1,lfp);
        lfp_hil_tf = hilbert(filtered_lfp);
        lfp_envelope = abs(lfp_hil_tf);
        smoothed_envelope = gaussfilt_2017(time,lfp_envelope,.004);
        env_zScore = zscore(smoothed_envelope);
        
        %ripple amplitude is the magnitude (zscore) of the smoothened envelope at
        %ripple peak
        midPt = round(length(lfp)/2);
        try
            spindleAmpResp(i).zScoreAmp(j,1) = max(env_zScore(midPt-2500:midPt+2500));
            
%         spindleAmpResp(i).zScoreAmp(j,1) = max(smoothed_envelope(midPt-150:midPt+150));
            spindleAmpResp(i).lfpTrace(j,:) = lfp;
%         spindleAmpResp(i).slowEnv(j,:) = env_zScore;
        
            %define time window (in imaging frames) surrounding ripple
            ripWindowIdx = spinLocs(j) - nFrames_window : spinLocs(j) + ...
                nFrames_window;

            spindleAmpResp(i).dff(j,:) = nanmean(dF(:,ripWindowIdx));
        end
end

%normalize the amplitude of ripples within each session
spindleAmpResp(i).zScoreAmp = spindleAmpResp(i).zScoreAmp./max(spindleAmpResp(i).zScoreAmp);
