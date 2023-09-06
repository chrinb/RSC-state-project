function plotEphys(sData)

% Plot ephys signals for visualization.

LFP_filtsig  = sData.ephysdata.lfpFilt;
signalLength = length(sData.ephysdata.lfp);
EMG_filtsig  = sData.ephysdata.EMGfilt;
sigma_LFP3   = sData.ephysdata.sigmaband;
sigma_LFP2   = (sData.ephysdata.sigmaband_2);
srate        = 2500;
delta_t      = 1/srate; 
tt           = delta_t:delta_t:signalLength*delta_t;


figure,
hAx(1) = subplot(611);
plot(tt, sData.daqdata.runSpeed)
set(gca, 'xlim', [0 600],'ylim', [0 10])
ylabel({'Running speed'})

hAx(2) = subplot(612);
plot(tt, EMG_filtsig)
set(gca, 'xlim', [0 600],'ylim', [-2 2])
ylabel({'EMG','0.1-100 KHz'})

hAx(3) = subplot(613);
plot(tt, sData.ephysdata3.lfp)
set(gca, 'xlim',[0 600])
ylabel({'S1 LFP'})

hAx(4) = subplot(614);
plot(tt, sigma_LFP3)
set(gca, 'xlim', [0 600],'ylim', [-.2 .2])
ylabel({'S1 LFP (sigma)','10-16 Hz'})

hAx(5) = subplot(615);
plot(tt, (sData.ephysdata2.lfp))
set(gca, 'xlim', [0 600])
ylabel({'RSC LFP'})

hAx(6) = subplot(616);
plot(tt, sigma_LFP2)
set(gca, 'xlim', [0 600],'ylim', [-.2 .2])
ylabel({'RSC LFP (sigma)','10-16 Hz'})

linkaxes(hAx);
sgt = sgtitle(sData.sessionInfo.sessionID);
sgt.FontSize = 8;