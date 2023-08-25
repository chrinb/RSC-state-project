function plotswrdff1(sData)

% Written by Christoffer Berge | Vervaeke Lab

% Plot mean DFF signal and all SWR time stamps (only use for awake recordings).

prompt = sprintf('All ripples? (y = yes | everything else = no) ');
allrip = input(prompt,'s');
frames = sData.daqdata.frame_onset_reference_frame;

if strcmp(allrip,'y')
    RippleIdx = sData.ephysdata.frameRipIdx; %keep all ripples
else
    prompt = sprintf('Only no movement ripples? (y = yes | everything else = no) ');
    riprun = input(prompt, 's');
    prompt = sprintf('Remove temporally close ripples? (y = yes | everything else = no) ');
    removerip = input(prompt, 's');
    if strcmp(riprun, 'y') && ~strcmp(removerip, 'y')
        RippleIdx = riprun2(sData); %removes all ripples if there is movement 2s before or after ripple peak
        RippleIdx = frames(RippleIdx);
    elseif strcmp(removerip, 'y') && ~strcmp(riprun, 'y')
        RippleIdx = RemoveRip(sData);
        RippleIdx = frames(RippleIdx);
    elseif strcmp(removerip, 'y') && strcmp(riprun, 'y')
        sData.ephysdata.frameRipIdx = riprun2(sData);
        RippleIdx = RemoveRip(sData);
        RippleIdx = frames(RippleIdx);
    end
end


sessionID = sData.sessionInfo.sessionID;

% find frames per second for imaging time vector
nrofFramesInRec = length(sData.imdata.roiSignals(2).newdff);
RecordingLengthSec = length(sData.ephysdata.lfp)/2500;
framesPerSec = nrofFramesInRec/RecordingLengthSec;

time = (0:length(sData.imdata.roiSignals(2).newdff)-1)/31;
time2 = (0:length(sData.ephysdata.lfp)-1)./2500;

AwakeSWR = zeros(1, length(sData.imdata.roiSignals(2).newdff(1,:)));
AwakeSWR(RippleIdx) = .1;


figure, 
subplot(211),
plot(time, mean(sData.imdata.roiSignals(2).newdff)), hold on, ...
plot(time,AwakeSWR),
xlabel('Time (s)'), 
ylabel('Mean DFF')
set(gca, 'xlim',[0 max(time)])
title(sessionID)

subplot(212),
plot(time2, sData.ephysdata.lfp)
set(gca, 'xlim', [0 max(time2)]);
ylabel('CA1 LFP (mV)')

