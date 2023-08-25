function rip_multi_swr_aligned_ecog(varargin)

% Calculate mean SWR-aligned RSC ECoG spectrogram. User inputs 

rsc_ecog  = varargin{1,1};

% Spectrogram parameters
window  = 128;              
overlap = 127;             
nFFT    = 50:10:300;       
Fs      = 2500;             



for n_ripples = 1:size(rsc_ecog,1)
    fprintf('Calculating spectrogram...ripple %d of %d\n',n_ripples,size(rsc_ecog,1));
    x = rsc_ecog(n_ripples,:);
    [~,~,~,cP] = spectrogram(x,window,overlap,nFFT,Fs);
    P(:,:,n_ripples) = cP;
end

P = mean(P,3);
zscoreP = zscore(P,0,2);
[~,F,T,~] = spectrogram(x,window,overlap,nFFT,Fs);
% Plot spectrogram
time = (-(length(T)/2):(length(T)/2)-1)./2500;
% spec_time = linspace(-1,1,length(T));

figure,
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