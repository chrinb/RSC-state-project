function plotswrecog(sData)

% Written by Christoffer Berge || Vervaeke lab

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
time = (0:length(sData.ephysdata.lfp)-1)/2500;
RippleIdx = sData.ephysdata.absRipIdx;

AwakeSWR = zeros(1, length(sData.ephysdata.lfp));
AwakeSWR(RippleIdx) = 1;

figure,
sgtitle(sessionID),
subplot(211)
plot(time, sData.ephysdata.lfp), hold on, ...
plot(time,AwakeSWR),
xlabel('Time (s)'), 
ylabel('mV')
if start > 0
    set(gca, 'xlim', [start stop], 'ylim', [-1 1]);
else
    set(gca, 'xlim',[0 max(time)], 'ylim', [-1 1])
end
title('CA1 LFP')

subplot(212)
plot(time, sData.ephysdata2.lfp), hold on, ...
plot(time,AwakeSWR),
xlabel('Time (s)'), 
ylabel('mV')
if start > 0
    set(gca, 'xlim', [start stop], 'ylim', [-.7 .7]);
else
    set(gca, 'xlim',[0 max(time)], 'ylim', [-.7 .7])
end
title('RSC ECoG')