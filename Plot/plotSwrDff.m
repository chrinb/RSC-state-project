function plotSwrDff(sData)

% Written by Christoffer Berge | Vervaeke Lab

% For sleep recordings with predefined sleep bouts, spindles and SWRs.
% Plots mean DF/F for bulk imaging data, NREM bouts, ECoG trace, and spindles.

signal = mean(sData.imdata.roiSignals(2).newdff);

prompt = sprintf('Select time interval? ');
signalSelect = input(prompt,'s');

if strcmp(signalSelect,'y')
    prompt = sprintf('Start (sec): ');
    start = input(prompt);
    prompt = sprintf('End (sec): ');
    stop = input(prompt);
else
    start = -1;
end
prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
spindle_band_select = input(prompt);

if spindle_band_select == 1
    ephys_select = [];
    spindle_select = [];
elseif spindle_band_select == 2
    ephys_select = num2str(2);
    spindle_select = '1016';
end
sessionID = sData.sessionInfo.sessionID;

temp_str = strcat('spindleStartEnd', num2str(spindle_select) );
temp_str2 = strcat('sigmaband', ephys_select);

ECoG_spindle = sData.ephysdata2.(temp_str)/2500;

[~, NREMstartend] = swrspindle(sData);
NREMstartend = NREMstartend./2500;

% find frames per second for imaging time vector
nrofFramesInRec = length(sData.imdata.roiSignals(2).newdff);
RecordingLengthSec = length(sData.ephysdata.lfp)/2500;
framesPerSec = nrofFramesInRec/RecordingLengthSec;

time_imaging = (0:length(sData.imdata.roiSignals(2).newdff)-1)/framesPerSec;
time_ephys = (0:length(sData.ephysdata2.lfp)-1)./2500;

%% Select SWRs for analysis
prompt = sprintf('All ripples? (y = yes | everything else = no) ');
allrip = input(prompt,'s');

frames = sData.daqdata.frame_onset_reference_frame;

if strcmp(allrip,'y') %keep all ripples
    unclassifiedSWRidx = sData.ephysdata.unclassified_swr;
    NREMspindleUncoupledSWRidx = sData.ephysdata.NREM_spindle_uncoupled_swr;
    NREMspindleCoupledSWRidx = sData.ephysdata.spindle_coupled_swr;
    
else
    prompt = sprintf('Remove locomotion SWR? (y = yes | everything else = no) ');
    riprun = input(prompt, 's');
    
    prompt = sprintf('Remove temporally close SWR? (y = yes | everything else = no) ');
    removerip = input(prompt, 's');
    
    % if remove locomotion SWR but not temporally close SWR
    if strcmp(riprun, 'y') && ~strcmp(removerip, 'y')
        [unclassifiedSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] = ripRunAn(sData); 
    % if remove temporally close SWR but not locomotion SWR
    elseif strcmp(removerip, 'y') && ~strcmp(riprun, 'y')
        [unclassifiedSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx] = removeCloseRip(sData);
    % if remove both temporally close and locotion SWR
    elseif strcmp(removerip, 'y') && strcmp(riprun, 'y')
        [unclassifiedSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] = ripRunAn(sData,1);
    end
end

AwakeSWR = zeros(1, length(sData.imdata.roiSignals(2).newdff(1,:)));
AwakeSWR(frames(unclassifiedSWRidx)) = .1;

UncoupledSWR = zeros(1, length(sData.imdata.roiSignals(2).newdff(1,:)));
UncoupledSWR(frames(NREMspindleUncoupledSWRidx)) = .1;

CoupledSWR = zeros(1, length(sData.imdata.roiSignals(2).newdff(1,:)));
CoupledSWR(frames(NREMspindleCoupledSWRidx)) = .1;

im_to_ephys_time = linspace(1, length(sData.ephysdata2.lfp), length(time_imaging));

figure,
sgtitle(sessionID),
subplot(211),
plot(time_imaging, signal), hold on, ...
plot(time_imaging,AwakeSWR),plot(time_imaging,UncoupledSWR),...
plot(time_imaging,CoupledSWR)
xlabel('Time (s)'), 
ylabel('Mean DFF')
set(gca, 'xlim',[0 max(time_imaging)])
set(gca, 'ylim',[min(signal)-.01, max(signal)+.01])
title('DF/F')


for i = 1:length(NREMstartend)
    x = [NREMstartend(i,1) NREMstartend(i,1) NREMstartend(i,2) NREMstartend(i,2)];
    y = [min(signal)-.01 max(signal)+.01 max(signal)+.01 min(signal)-.01];
    patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .3);
end
if start > 0
    set(gca, 'xlim', [start stop]);
else
end

subplot(212)
plot(time_ephys, sData.ephysdata2.lfp), hold on,
plot(time_ephys, (sData.ephysdata2.(temp_str2)*3)+.1)
xlabel('Time (s)'), 
ylabel('mV')
set(gca, 'xlim',[0 max(time_ephys)])
set(gca, 'ylim', [-1,1])
title('RSC ECoG')

for i = 1:length(NREMstartend)
    x = [NREMstartend(i,1) NREMstartend(i,1) NREMstartend(i,2) NREMstartend(i,2)];
    y = [-2 2 2 -2];
    patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .3);
end
if start > 0
    set(gca, 'xlim', [start stop]);
else
end

for k = 1:length(sData.ephysdata2.(temp_str))
    z = [ECoG_spindle(k,1) ECoG_spindle(k,1) ECoG_spindle(k,2) ECoG_spindle(k,2) ];
    v = [-2 2 2 -2];
    patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .5);
end

if start > 0
    set(gca, 'xlim', [start stop]);
else
end