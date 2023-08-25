function r_swr_emg = event_cross_corr_awake(varargin)

% Calculate mean cross-correlation between SWR, spindle, and delta power
% time series. 

sData                   = varargin{1,1};
LFP                     = sData.ephysdata.lfp;
ECoG                    = sData.ephysdata2.lfp;
% sigma_hilbert           = hilbert(sData.ephysdata2.sigmaband);
% sigma_pow               = abs(sigma_hilbert);
swr_hilbert             = hilbert(sData.ephysdata.ripplefreq);
swr_pow                 = abs(swr_hilbert);
srate                   = 2500;
% single_swr_idx          = sData.ephysdata.spindle_coupled_swr;  
% spindle_coupled_swr_idx = sData.ephysdata.NREM_spindle_uncoupled_swr;  
swr_idx                 = sData.ephysdata.absRipIdx;
time_vec                = (-4:1/srate:4);
% NREM_spindle_idx        = sData.ephysdata2.NREMAbsSpindleIdx;
window                  = 4;
emg_filt                = sData.ephysdata3.EMGfilt;
emg_pow_env             = abs( hilbert(emg_filt));
run_speed               = sData.daqdata.runSpeed;
% [nrem_SO, nrem_delta_waves] = mark_slow_wave(sData);
% nrem_delta_waves = slow_wave_an2(sData);
delta_hilbert = hilbert(sData.ephysdata2.deltaband);
% delta_pow     = abs(delta_hilbert);
% delta_idx     = vertcat(nrem_delta_waves(:,3), nrem_SO(:,3));

%% SWR-spindle & SWR-delta cross correlation

% Preallocate
r_swr_spindle = zeros(length(swr_idx),srate*window*window+1);
r_swr_delta   = zeros(length(swr_idx),srate*window*window+1);
r_swr_emg     = zeros(length(swr_idx),srate*window*window+1);
r_swr_speed   = zeros(length(swr_idx),srate*window*window+1);
% Loop over nr of SWRs
for swr_nr = 1:length(swr_idx)
    time_win_start = swr_idx(swr_nr) - window*srate;
    time_win_end   = swr_idx(swr_nr) + window*srate;
    
    if time_win_start > 1 && time_win_end < length(LFP)
        swr_pow_win     = swr_pow(time_win_start:time_win_end);
        swr_pow_win     = swr_pow_win - mean(swr_pow_win);
        
%         spindle_pow_win = sigma_pow(time_win_start:time_win_end);
%         spindle_pow_win = spindle_pow_win - mean(spindle_pow_win);
%         
%         delta_pow_win   = delta_pow(time_win_start:time_win_end);
%         delta_pow_win   = delta_pow_win - mean(delta_pow_win);

        emg_pow_win     = emg_pow_env(time_win_start:time_win_end);
        emg_pow_win     = emg_pow_win - mean(emg_pow_win);

        run_speed_win   = run_speed(time_win_start:time_win_end);
        run_speed_win   = run_speed_win - mean(run_speed_win);

%         r_swr_spindle(swr_nr,:) = xcorr(swr_pow_win, spindle_pow_win, 'normalized');
%         r_swr_delta(swr_nr,:)   = xcorr(swr_pow_win, delta_pow_win, 'normalized');
        r_swr_emg(swr_nr,:)     = xcorr(swr_pow_win,emg_pow_win,'normalized');
        r_swr_speed(swr_nr,:)   = xcorr(swr_pow_win, run_speed_win,'normalized');
    end
end

%% Delta-spindle & delta-SWR cross-correlation
% r_delta_spindle = zeros(length(delta_idx),srate*window*window+1);
% r_delta_swr     = zeros(length(delta_idx),srate*window*window+1);
% 
% for delta_nr = 1:length(delta_idx)
%     time_win_start = delta_idx(delta_nr) - window*srate;
%     time_win_end   = delta_idx(delta_nr) + window*srate;
%     
%     if time_win_start > 1 && time_win_end < length(LFP)
%         spindle_pow_win2 = sigma_pow(time_win_start:time_win_end);  
%         spindle_pow_win2 = spindle_pow_win2 - mean(spindle_pow_win2);
%         
%         delta_pow_win2  = delta_pow(time_win_start:time_win_end);
%         delta_pow_win2  = delta_pow_win2 - mean(delta_pow_win2);
%         
%         swr_pow_win2    = swr_pow(time_win_start:time_win_end);
%         swr_pow_win2    = swr_pow_win2 - mean(swr_pow_win2);
% 
%         r_delta_spindle(delta_nr,:) = xcorr(delta_pow_win2, spindle_pow_win2,'normalized');
%         r_delta_swr(delta_nr,:) = xcorr(delta_pow_win2, swr_pow_win2,'normalized');
% 
%     end
% end
% 
% %% Spindle-SWR cross correlation
% r_spin_swr = zeros(length(NREM_spindle_idx),srate*window*window+1);
% 
% % Loop over nr of NREM spindles
% for spin_nr = 1:length(NREM_spindle_idx)
%     time_win_start = round( NREM_spindle_idx(spin_nr) - window*srate);
%     time_win_end   = round( NREM_spindle_idx(spin_nr) + window*srate);
%     
%     if time_win_start > 1 && time_win_end < length(LFP)
%         swr_pow_win2     = swr_pow(time_win_start:time_win_end);
%         swr_pow_win2     = swr_pow_win2 - mean(swr_pow_win2);
%         
%         spindle_pow_win2 = sigma_pow(time_win_start:time_win_end);
%         spindle_pow_win2 = spindle_pow_win2 - mean(spindle_pow_win2);
%         
%         delta_pow_win2 = delta_pow(time_win_start:time_win_end);
%         delta_pow_win2 = delta_pow_win2 - mean(delta_pow_win2);
%         
%         r_spin_swr(spin_nr,:) = xcorr(spindle_pow_win2, swr_pow_win2,'normalized');
%         r_spin_delta(spin_nr,:) = xcorr(spindle_pow_win2, delta_pow_win2, 'normalized');
%     end
% end
% 
% mean_swr_spindle_xcorr   = mean(r_swr_spindle);
% mean_swr_delta_xcorr     = mean(r_swr_delta);
% mean_delta_spindle_xcorr = mean(r_delta_spindle);
% mean_delta_swr_xcorr     = mean(r_delta_swr);
% mean_spindle_swr_xcorr   = mean(r_spin_swr);
% mean_spindle_delta_xcorr = mean(r_spin_delta);
% mean_swr_emg_xcorr       = mean(r_swr_emg);

% mean_xcorr = {mean_swr_spindle_xcorr, mean_swr_delta_xcorr,mean_delta_spindle_xcorr,...
%     mean_delta_swr_xcorr, mean_spindle_swr_xcorr, r_swr_emg};
% 
% %% Stuff for plotting
% time1                  = -8:1/2500:8;
% SE_spindle_swr_xcorr   = std(r_spin_swr) ./ sqrt(size(r_spin_swr,1)); 
% SE_swr_delta_xcorr     = std(r_swr_delta) ./ sqrt(size(r_swr_delta,1)); 
% SE_spindle_swr_xcorr   = std(r_spin_swr) ./ sqrt(size(r_spin_swr,1)); 
% SE_delta_spindle_xcorr = std(r_delta_spindle) ./ sqrt(size(r_delta_spindle,1)); 

