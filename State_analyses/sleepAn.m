function sleepAn(sData)

EMG = sData.ephysdata.lfp;
LFP2 = sData.ephysdata2.lfp;
LFP3 = sData.ephysdata3.lfp;
signalLength = length(sData.ephysdata.lfp);
shape = [0 0 1 1 0 0];
transw = 0.1;
srate = 2500;
nyquist = srate/2;
delta_t = 1/srate; 
tt = delta_t:delta_t:signalLength*delta_t;


LFP3_filtsig = sData.ephysdata.lfp3Filt;
LFP2_filtsig = sData.ephysdata.lfp2Filt;
EMG_filtsig = sData.ephysdata.EMGfilt;
delta_filtsig = sData.ephysdata.deltaband2;
theta_filtsig = sData.ephysdata.thetaband2;
sigma_filtsig = sData.ephysdata.sigmaband2;
delta_2_filtsig = (sData.ephysdata.deltaband3);
theta_2_filtsig = (sData.ephysdata.thetaband3);
sigma_2_filtsig = (sData.ephysdata.sigmaband3);
%ripple_filtsig = sData.ephysdata.ripplefreq;

hilbert_tf_EMG = hilbert(EMG_filtsig);
%hilbert_tf_delta = hilbert(delta_filtsig);
%hilbert_tf_theta = hilbert(theta_filtsig);
%hilbert_tf_sigma = hilbert(sigma_filtsig);

k = 1;
for i = 1:(1*srate):(signalLength-(5*srate))
    delta_mean(k) = rms(delta_filtsig(i:(i+(5*srate))));
    theta_mean(k) = rms(theta_filtsig(i:(i+(5*srate))));
    sigma_mean(k)= rms(sigma_filtsig(i:(i+(5*srate))));
    ratiotheta(k) = theta_mean(k)/(theta_mean(k)+delta_mean(k));
    
    time_vec(k) = tt(i);
    k = k+1;
end
time_vec2 = time_vec';

p = 1;
for i = 1:(1*srate):(signalLength-(5*srate))
    delta_2_mean(p) = rms(delta_2_filtsig(i:(i+(5*srate))));
    theta_2_mean(p) = rms(theta_2_filtsig(i:(i+(5*srate))));
    sigma_2_mean(p)= rms(sigma_2_filtsig(i:(i+(5*srate))));
    ratiotheta2(p) = theta_2_mean(p)/(theta_2_mean(p)+delta_2_mean(p));
    
    time_vec(p) = tt(i);
    p = p+1;
end
time_vec2 = time_vec';

figure,
subplot(711)
plot(tt, sData.daqdata.runSpeed)
set(gca, 'xlim', [0 1800],'ylim', [0 10])
ylabel({'Running speed'})

subplot(712)
plot(tt, LFP3_filtsig)
set(gca, 'xlim', [0 1800],'ylim', [-.5 .5])
ylabel({'LFP (3)','0.5-30 Hz'})

subplot(713)
plot(tt, EMG_filtsig)
set(gca, 'xlim', [0 1800],'ylim', [-2 2])
ylabel({'EMG','0.1-100 KHz'})

subplot(714)
plot(delta_mean), set(gca, 'ylim', [0, 0.1])
ylabel({'Delta','0.5-4 Hz'})

subplot(715)
plot(theta_mean)
ylabel({'Theta','5-9 Hz'})

subplot(716)
plot(sigma_mean)
ylabel({'Sigma','10-16 Hz'})

subplot(717)
plot(ratiotheta)
ylabel({'Theta/(Theta+Delta)'}), yline(0.5, '--');

sgt = sgtitle(sData.sessionInfo.sessionID);
sgt.FontSize = 8;

figure,
subplot(711)
plot(tt, sData.daqdata.runSpeed)
set(gca, 'xlim', [0 1800],'ylim', [0 10])
ylabel({'Running speed'})

subplot(712)
plot(tt, LFP2_filtsig)
set(gca, 'xlim', [0 1800], 'ylim', [-.5 .5])
ylabel({'LFP (2)','0.5-30 Hz'})

subplot(713)
plot(tt, EMG_filtsig)
set(gca, 'xlim', [0 1800],'ylim', [-2 2])
ylabel({'EMG','0.1-100 KHz'})

subplot(714)
plot(delta_2_mean), set(gca, 'ylim', [0, 0.1])
ylabel({'Delta','0.5-4 Hz'})

subplot(715)
plot(theta_2_mean)
ylabel({'Theta','5-9 Hz'})

subplot(716)
plot(sigma_2_mean)
ylabel({'Sigma','10-16 Hz'})

subplot(717)
plot(ratiotheta2)
ylabel({'Theta/(Theta+Delta)'}), yline(0.5, '--');


sgt = sgtitle(sData.sessionInfo.sessionID);
sgt.FontSize = 8;