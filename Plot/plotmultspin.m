function plotmultspin(sData, x)

prompt = sprintf('Choose LFP signal ( 2 = LFP2, 3 = LFP3, 4 = LFP4): ');
signal = input(prompt);

if signal == 2
        aa = sData.ephysdata2.absSpindleIdx;

elseif signal == 3
        aa = sData.ephysdata3.absSpindleIdx;

elseif signal == 4
        aa = sData.ephysdata4.absSpindleIdx;
end

% compute spindle window

spinWinStart = round( aa(x) - (2500*4) );
spinWinEnd = round( aa(x) + (2500*4) );

time = (0:20001-1)/ 2500;
window  = 2500;              
overlap = 2400;              
nFFT    = 5:.1:20;         
fs = 2500;


figure, 
subplot(611),
plot(time, sData.ephysdata2.lfp(spinWinStart:spinWinEnd))
set(gca, 'ylim',[-.5 .5])
title('RSC')

x = sData.ephysdata2.lfp(spinWinStart:spinWinEnd);
[s,F,T,P] = spectrogram(x,window,overlap,nFFT,fs);

subplot(612)
contourf(T,F, (abs(P)),200,'edgecolor','none');
colormap(jet)
xlabel('Time (s)')

subplot(613),plot(time, sData.ephysdata3.lfp(spinWinStart:spinWinEnd))
set(gca, 'ylim',[-.5 .5])
title('PFC')

x = sData.ephysdata3.lfp(spinWinStart:spinWinEnd);
[s,F,T,P2] = spectrogram(x,window,overlap,nFFT,fs);

subplot(614)
contourf(T,F, (abs(P2)),200,'edgecolor','none');
colormap(jet)
xlabel('Time (s)')


subplot(615),plot(time, sData.ephysdata4.lfp(spinWinStart:spinWinEnd))
set(gca, 'ylim',[-.5 .5])
title('S1')

x = sData.ephysdata4.lfp(spinWinStart:spinWinEnd);
[s,F,T,P3] = spectrogram(x,window,overlap,nFFT,fs);

subplot(616)
contourf(T,F, (abs(P3)),200,'edgecolor','none');
colormap(jet)
xlabel('Time (s)')



