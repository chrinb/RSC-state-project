function [] = plotSWRECOGspectrogram(sData)

sessionID = sData.sessionInfo.sessionID;

window  = 128;              % Window size for computing the spectrogram (FFT) [# samples]
overlap = 127;              % Overlap of the windows for computing the spectrogram [# samples]
nFFT    = 50:10:300;        % Vector defining the frequencies for computing the FFT
Fs      = 2500;             % Signal sampling frequency.

rippleIdx = sData.ephysdata.absRipIdx;
nRipples = length(sData.ephysdata.absRipIdx);
ECoG = sData.ephysdata2.lfp;

% find RSC ECoG snippets
for k = 1:nRipples
    ECoG_snippet_start = rippleIdx(k)-750;
    ECoG_snippet_end = rippleIdx(k)+750;
    ECoG_snippet(k,:) = ECoG(ECoG_snippet_start:ECoG_snippet_end);
end

for z = 1:size(ECoG_snippet,1)
    fprintf('Calculating spectrogram...ripple %d of %d\n',z,size(ECoG_snippet,1));
    x = ECoG_snippet(z,:);
    [~,~,~,cP] = spectrogram(x,window,overlap,nFFT,Fs);
    P(:,:,z) = cP;
end
P = mean(P,3);
zscoreP = zscore(P,0,2);
[~,F,T,~] = spectrogram(x,window,overlap,nFFT,Fs);
% Plot spectrogram
time = (-(length(T)/2):(length(T)/2)-1)./2500;

figure,
sgtitle(sessionID),
subplot(211)
contourf(time,F,10*log10(abs(P)),200,'edgecolor','none');
colormap(jet)
colorbar
title('Avg HPC-SWR aligned ECoG spectrogram')
ylabel('Frequency (Hz)')

subplot(212)
contourf(time,F,10*log10(abs(zscoreP)),200,'edgecolor','none');
colormap(jet)
caxis([-5 5])
colorbar
title('Avg HPC-SWR aligned ECoG spectrogram (z-score)')
ylabel('Frequency (Hz)')
xlabel('Time from HPC SWR peak (s)')