function plot_bulk_spindle_dff(sData)

%  Written by Christoffer Berge | Vervaeke Lab

% Plot mean DF/F signal for bulk imaging sessions, NREM bouts, and sleep
% spindles. 

% Get session name
sessionID = sData.sessionInfo.sessionID;

% Find start/end of NREM bouts and convert to seconds.
[~, NREMstartend] = swrspindle(sData);
NREMstartend      = NREMstartend./2500;

% find frames per second for imaging time vector
nrofFramesInRec    = length(sData.imdata.roiSignals(2).newdff);
RecordingLengthSec = length(sData.ephysdata.lfp)/2500;
framesPerSec       = nrofFramesInRec/RecordingLengthSec;

% Create time vector (in seconds) for plotting
time = (0:length(sData.imdata.roiSignals(2).newdff)-1)/framesPerSec;
time2 = (0:length(sData.ephysdata2.lfp)-1)./2500;

% Find start/end of sleep spindles (in seconds)
ECoG_spindle = sData.ephysdata2.NREMspindleStartEnd/2500;

% Mark center time point of sleep spindles
spindle_center                                   = zeros(1, length(sData.imdata.roiSignals(2).newdff(1,:)));
spindle_temp                                     = sData.ephysdata2.absSpindleIdx./2500;
spindle_center(round(spindle_temp*framesPerSec)) = .1;

% Test to see random non-spindle NREM windows
[~,rand_window] = rand_nrem_times(sData);

%% Plot results
figure
subplot(211)
plot(time, mean(sData.imdata.roiSignals(2).newdff)), hold on, ...
plot(time,spindle_center)
xlabel('Time (s)'), 
ylabel('Mean DFF')
set(gca, 'xlim',[0 max(time)], 'ylim',[-.4 .4])
title(sessionID)

% Loop over NREM bouts and plot them as individual patches
for i = 1:length(NREMstartend)
    x = [NREMstartend(i,1) NREMstartend(i,1) NREMstartend(i,2) NREMstartend(i,2)];
    y = [-.6 .6 .6 -.6];
    h = patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .3);
end

% Loop over sleep spindles and plot them as individual patches
for k = 1:length(sData.ephysdata2.NREMspindleStartEnd)
    z = [ECoG_spindle(k,1) ECoG_spindle(k,1) ECoG_spindle(k,2) ECoG_spindle(k,2) ];
    v = [-.6 .6 .6 -.6];
    h2 = patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .5);
end

for o = 1:length(rand_window)
    r = [rand_window(o,1) rand_window(o,1) rand_window(o,2) rand_window(o,2) ]./2500;
    e = [-.6 .6 .6 -.6];
    h3 = patch(r,e , 'red', 'edgecolor', 'none', 'FaceAlpha', .5);
end
legend( [h, h2, h3], 'NREM','Sleep spindles', 'random');

subplot(212)
plot(time2, sData.daqdata.runSpeed)
ylabel('Running speed cm/s')
xlabel('Time (s)')
set(gca, 'xlim',[0 time2(end)])