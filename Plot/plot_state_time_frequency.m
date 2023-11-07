function [freq_ecog, time_ecog, pow_ecog, freq_lfp, time_lfp, pow_lfp] = plot_state_time_frequency(sData)

% Written by Christoffer Berge | Veraeke lab

% Function that computes time-frequency data from ECoG channel

%% Define variables
ECoG        = downsample(sData.ephysdata2.lfp, 10);
LFP         = downsample(sData.ephysdata.lfp, 10);
% ECoG        = sData.ephysdata2.lfp;
% LFP         = sData.ephysdata.lfp;
srate       = 250;
% time_vector = (1:length(ECoG))/srate;
window      = 300;              
overlap     = 280;  
nFFT        = 1:.5:30;        

% EcoG
[~ ,freq_ecog, time_ecog, pow_ecog] = spectrogram(ECoG,window,overlap,nFFT,srate);

% HPC CA1 LFP
[~ ,freq_lfp, time_lfp, pow_lfp] = spectrogram(LFP,window,overlap,nFFT,srate);

% figure
% clf
% figure(1), contourf(time_ecog, freq_ecog, 10*log10( pow_ecog),300,'linecolor','none'), 
% clim([-45 -20]), colormap jet

% figure(2), imagesc(time_ecog,flipud(freq_ecog), 10*log10(pow_ecog) ), colorbar
% clim([-45 -20]), colormap jet




% clf, contourf(time_instants,cFrequencies,10*log10( pow_spectrum),300,'linecolor','none'), colorbar, caxis([-50 -20]), colormap jet
% xlabel('Time (s)'), ylabel('Frequency (Hz)')
% figure(1), clf,  imagesc(time_vector, cFrequencies, abs(pow_spectrum) )