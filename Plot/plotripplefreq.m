function plotripplefreq(sData)

% Written by Christoffer Berge | Vervaeke Lab

swrIdx = zeros(1, length(sData.ephysdata.lfp));
swrIdx(sData.ephysdata.absRipIdx) = 0.02;

time = (0:length(sData.ephysdata.lfp)-1)/2500;

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

sessionID = sData.sessionInfo.sessionID;
LFP = sData.ephysdata.lfp;
ECoG = sData.ephysdata2.lfp;
CA1_ripple_power = (2*abs(sData.ephysdata.ripplefreq)).^2;
RSC_ripple_power = (2*abs(sData.ephysdata2.ripplefreq)).^2;

swrIdx = zeros(1, length(sData.ephysdata.lfp));
swrIdx2 = zeros(1, length(sData.ephysdata.lfp));
swrIdx3 = zeros(1, length(sData.ephysdata.lfp));
swrIdx4 = zeros(1, length(sData.ephysdata.lfp));
swrIdx(sData.ephysdata.absRipIdx) = 0.04;
swrIdx2(sData.ephysdata.absRipIdx) = 0.02;
swrIdx3(sData.ephysdata.absRipIdx) = 0.3;
swrIdx4(sData.ephysdata.absRipIdx) = 0.;

time = (0:length(sData.ephysdata.lfp)-1)/2500;

figure,
hAx(1) = subplot(211);
plot(time, LFP*.2)
title(sessionID)
set(gca, 'xlim', [0 max(time)]);
ylabel('CA1 LFP (mV)')
if start > 0
    set(gca, 'xlim', [start stop]);
else
end

hold on
plot(time, CA1_ripple_power -.3), hold on, plot(time, swrIdx -.3, 'color', 'k','linewidth',1)
set(gca, 'xlim', [0 max(time)]);
ylabel('Power (mV^2)')
if start > 0
    set(gca, 'xlim', [start stop]);
else
end

hAx(2) = subplot(212);
plot(time, ECoG*.2)
set(gca, 'xlim', [0 max(time)], 'ylim',[-.5 .5]);
ylabel('RSC ECoG (mV)')
if start > 0
    set(gca, 'xlim', [start stop]);
else
end

hold on
plot(time, RSC_ripple_power*10-.1), hold on, plot(time, swrIdx-.05, 'linewidth',1)
set(gca, 'xlim', [0 max(time)], 'ylim',[-.1 .1]);

if start > 0
    set(gca, 'xlim', [start stop]);
else
end
ylabel('Power (mV^2)')
xlabel('Time (s)')

linkaxes(hAx, 'x');