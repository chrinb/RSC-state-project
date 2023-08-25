function plotswrspindle(sData, varargin)

% Written by Christoffer Berge || Vervaeke lab

% plots raw HPC CA1 LFP signal with detected SWRs (classified as awake, NREM
% spindle-uncoupled SWR and spindle-coupled SWR), RSC ECoG signal with 
% detected spindles. 
sessionID = sData.sessionInfo.sessionID;
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


LFP = sData.ephysdata.lfp;
ECoG = sData.ephysdata2.lfp;
SWR_freq_power = abs(hilbert(sData.ephysdata.ripplefreq).^2);
spindle_freq_power = abs(hilbert(sData.ephysdata2.sigmaband).^2);

awake_swr = sData.ephysdata.unclassified_swr;
NREM_spindle_uncoupled_swr = sData.ephysdata.NREM_spindle_uncoupled_swr;
spindle_coupled_swr = sData.ephysdata.spindle_coupled_swr;
spindle_idx = round(sData.ephysdata2.absSpindleIdx);

[~, NREMstartend] = swrspindle(sData);

% find spindle times
spindle = zeros(1, length(sData.ephysdata2.lfp));
spindle(round(sData.ephysdata2.absSpindleIdx)) = .1;

% find SWR times & spindle centre 
allSWR = zeros(1, length(sData.ephysdata.lfp));
allSWR(sData.ephysdata.absRipIdx) = .02;

awakeSWR = zeros(1, length(sData.ephysdata.lfp));
awakeSWR(awake_swr) = 1;

NREM_spindle_uncoupledSWR = zeros(1, length(sData.ephysdata.lfp));
NREM_spindle_uncoupledSWR(NREM_spindle_uncoupled_swr) = 1;

spindle_coupledSWR = zeros(1, length(sData.ephysdata.lfp));
spindle_coupledSWR(spindle_coupled_swr) = 1;

spindleIdx = zeros(1, length(sData.ephysdata.lfp));
spindleIdx(spindle_idx) = 1;

NREMstartend = NREMstartend./2500;
time = (0:length(sData.ephysdata.lfp)-1) / 2500;

ECoG_spindle = sData.ephysdata2.spindleStartEnd/2500;

figure, 
hAx(1) = subplot(511);
sgtitle(sessionID) 
plot(time, LFP), set(gca,'xlim', [0 max(time)]), hold on,
plot(time, awakeSWR, 'linewidth', 1), 
plot(time, NREM_spindle_uncoupledSWR, 'linewidth', 1),
plot(time, spindle_coupledSWR, 'linewidth', 1)
set(gca, 'ylim',[-1 1])
xlabel('Time (s)')
ylabel('mV')
title('HPC LFP')
%legend('Raw data','awake swr', 'NREM spindle-uncoupled SWR', 'spindle-coupled SWR')

for i = 1:length(NREMstartend)
    x = [NREMstartend(i,1) NREMstartend(i,1) NREMstartend(i,2) NREMstartend(i,2)];
    y = [-1.5 1.5 1.5 -1.5];
    patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .3);
end
if start > 0
    set(gca, 'xlim', [start stop]);
else
end

hAx(2) = subplot(512);
plot(time, SWR_freq_power,time,allSWR), set(gca,'xlim', [0 max(time)]), 
xlabel('Time (s)')
ylabel('Power (mV^2)')
set(gca, 'ylim',[0 0.04])
if start > 0
    set(gca, 'xlim', [start stop]);
else
end

hAx(3) = subplot(513);
plot(time, ECoG, time,spindle),  set(gca, 'ylim',[-.6 .5], 'xlim', [0 max(time)]), 
xlabel('Time (s)')
ylabel('mV')
title('RSC')
hold on

for i = 1:length(NREMstartend)
    x = [NREMstartend(i,1) NREMstartend(i,1) NREMstartend(i,2) NREMstartend(i,2)];
    y = [-.6 .6 .6 -.6];
    patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .3);
end

    for k = 1:length(sData.ephysdata2.spindleStartEnd)
    z = [ECoG_spindle(k,1) ECoG_spindle(k,1) ECoG_spindle(k,2) ECoG_spindle(k,2) ];
    v = [-.6 .6 .6 -.6];
    patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .5);
    end
if start > 0
    set(gca, 'xlim', [start stop]);
else
end

hAx(4) = subplot(514);
plot(time, spindle_freq_power), set(gca,'xlim', [0 max(time)]), hold on,
plot(time, spindleIdx, 'linew',1)
xlabel('Time (s)')
ylabel('Power (mV^2)')
set(gca, 'ylim',[0 0.04])
if start > 0
    set(gca, 'xlim', [start stop]);
else
end


hAx(5) = subplot(515);
plot(time, sData.ephysdata2.lfpFilt),hold on
plot(time, sData.ephysdata2.deltaband, 'linew',1)
set(gca, 'xlim', [time(1), time(end)])
linkaxes(hAx,'x');

