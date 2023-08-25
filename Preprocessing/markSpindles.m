function [sData] = markSpindles(sData)


srate = 2500;
%freqfilter = [10 16];
lfpsignal = sData.ephysdata.lfp;
%lfpsignal2 = sData.ephysdata3.lfp;

hz = linspace(0, srate, floor(length(srate)/2+1));

nSnips = floor(length(lfpsignal)/(srate)) - 1;
time = linspace(0,length(lfpsignal),length(lfpsignal))/(srate);
timeRound = round(time,3);
spindleLocs = [];
runSignal = sData.daqdata.runSpeed;

window_size = 1;
window_size_index = window_size * srate;

%freqL = freqFilter(1);
%freqU = freqFilter(2);
nyquistFs = srate/2;

U_threshold = 1.5;

LFP2sigma = sData.ephysdata2.sigmaband2;
% LFP3sigma = sData.ephysdata.sigmaband3;
% LFP4sigma = sData.ephysdata.sigmaband4;

% Get sigma frequency band amplitude envelope
LFP2__hilbert = hilbert(LFP2sigma);
LFP2_envelope = abs(LFP2__hilbert);

% LFP3_hilbert = hilbert(LFP3sigma);
% LFP4_envelope = abs(LFP3_hilbert);

% LFP4_hilbert = hilbert(LFP4sigma);
% LFP4_envelope = abs(LFP4_hilbert);

% Index NREM episodes

%thresholdSigma

% Calculate upper and lower threshold of the sigma envelope from NREM
% episodes
LFP2LowThres = mean(LFP2_envelope) + 0.2*std(LFP2_envelope);

[~,locs,~,~] = findpeaks(LFP2smoothed_envelop-upper_thresh,srate,'MinPeakDistance',0.025,'WidthReference','halfhprom','Annotate','extents','WidthReference','halfprom');
spindleLocs = round(locs,3); 

spindleSnips = struct();
spindleIdx = zeros(1,length(spindleLocs));
spindleLocs = round(spindleLocs,3);

%convert the ripple locations from time to sample
for i = 1:length(spindleLocs)
        lfpPeakIdx = find(timeRound == spindleLocs(i));
        lfpStartIdx = lfpPeakIdx(1) - (0.5*srate);
        %if the ripple timepoint is near the beginning of the trace
        if lfpStartIdx < 0; lfpStartIdx = 1; end
        lfpEndIdx = lfpPeakIdx(1) + (0.5*srate);
        %if the ripple timepoint is near the end of the trace
        if lfpEndIdx > length(lfpsignal); lfpEndIdx = length(lfpsignal); end
        spindleSnips(i).lfp = lfpsignal(lfpStartIdx:lfpEndIdx);
        spindleIdx(i) = lfpPeakIdx(1);
    if runSignal(spindleIdx(i)) ~= 0; spindleIdx(i) = NaN; end %take out timepoints when animal is walking
%     if runSignal(moving_mean_move(i)) ~= 0; rippleIdx(i) = NaN; end %take out timepoints when animal is walking
    
end

%remove NaNs and timepoints too close together (likely identifying the same ripple
%waveform)
spindleSnips(isnan(spindleIdx)) = [];
spindleIdx(isnan(spindleIdx)) = [];



[final_rippleLFP,final_spindleLocs] = inspectSpindles(spindleSnips,spindleIdx,lfpsignal);
sData.ephysdata.absSpindleIdx = final_spindleLocs;
sData.ephysdata.spindleSnips = final_rippleLFP;
try frames = sData.daqdata.frame_onset_reference_frame;
catch frames = sData.daqdata.frameIndex;
end
sData.ephysdata.frameSpindleIdx = frames(sData.ephysdata.absSpindleIdx);
sData.ephysdata.SpindleEnvThr = U_threshold;