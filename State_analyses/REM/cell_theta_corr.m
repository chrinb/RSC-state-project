function [sorted_corr_ampl, sorted_corr_idx_ampl] = cell_theta_corr(sData, params)

% Written by Christoffer Berge | Vervaeke lab

% Calculate the Pearson correlation between cells and theta rhythm in h
% hippocampal LFP and neocortical ECoG signal during REM sleep.

% The code tries a few different strategies: Upsampling DF/F traces to
% match the sample rate of ephys (i.e. 2500Hz), or downsampling ephys
% signal to 2P framerate (~31Hz). 

% Also, try correlate theta band signal and theta band power amplitude
% envelope.


% Select cell type to analyze
switch params.cell_type

    case 'in'
        [~, in] = remove_cells_longitudinal(sData);
        dff     = sData.imdata.roiSignals(2).newdff(in, :);

    case 'pc'
        [pc, ~] = remove_cells_longitudinal(sData);
        dff     = sData.imdata.roiSignals(2).newdff(pc, :);

    case 'axon'
        dff = sData.imdata.roiSignals(2).mergedAxonsDffFilt;
end

%% Get raw and theta band filtered LFP and ECoG, sample rates, and
% ephys-to-imaging variable
lfp                 = sData.ephysdata.lfp;
ecog                = sData.ephysdata2.lfp;
lfp_theta           = sData.ephysdata.thetaband;
ecog_theta          = sData.ephysdata2.thetaband;
ephys_srate         = 2500;
imaging_srate       = find_imaging_framerate(sData);
frames              = sData.daqdata.frame_onset_reference_frame;

% Get REM sleep start/stop times
rem_theta_times = rem_sleep(sData);

% Preallocate
[time_ephys, time_imaging, dff_rem, rem_lfp, rem_ecog, rem_theta_lfp, rem_theta_ecog,...
    rem_theta_ampl, rem_ecog_ampl]  = deal( cell( size(rem_theta_times,1)));

% Loop over REM episodes in recording and extract snippets of raw and
% theta filtered data, and corresponding DF/F traces 
for rem_nr = 1:size(rem_theta_times,1)
    
    rem_snippet_ephys     = rem_theta_times(rem_nr, 1):rem_theta_times(rem_nr, 2);
    rem_snippet_imaging   = frames(rem_theta_times(rem_nr, 1)):frames(rem_theta_times(rem_nr, 2));

    time_ephys{rem_nr}    = (0:length(rem_snippet_ephys)-1)./ephys_srate;
    time_imaging{rem_nr}  = (0:length(rem_snippet_imaging)-1)./imaging_srate;

    rem_lfp{rem_nr}        = lfp(rem_snippet_ephys);
    rem_theta_lfp{rem_nr}  = lfp_theta(rem_snippet_ephys);
    rem_theta_ampl{rem_nr}   = abs( hilbert(lfp_theta(rem_snippet_ephys)));

    rem_ecog{rem_nr}       = ecog(rem_snippet_ephys);
    rem_theta_ecog{rem_nr} = ecog_theta(rem_snippet_ephys);

    dff_rem{rem_nr,:}      = dff(:, rem_snippet_imaging);
end

%% Strategy 1: Downsample ephys signal to match imaging data

% Downsample ephys data to imaging sample rate
downsample_factor = ephys_srate/imaging_srate;

% Lowpass filter signal at new Nyquist frequency
fkern        = fir1( round(14*imaging_srate/2), (imaging_srate/2) / (ephys_srate/2) );
f_theta      = filtfilt( fkern, 1, rem_theta_lfp{1});
f_theta_ampl = filtfilt( fkern, 1, rem_theta_ampl{1});

% Downsampled theta
ds_theta = f_theta(1:downsample_factor:end);

ds_theta_ampl = f_theta_ampl(1:downsample_factor:end);
% newtime = (0:length(ds_theta)-1)/ imaging_srate;

% Preallocate correlation coefficients
[theta_corr_ds, theta_ampl_corr_ds] = deal( zeros( size(dff,1), 1));

for roi_nr = 1:size(dff_rem{1}, 1)
    
    % Pearson correlation between downsampled theta and DF/F for each ROI
    tmp_corr_mat             = corrcoef(ds_theta', dff_rem{1}(roi_nr, :));
    theta_corr_ds(roi_nr, 1) = tmp_corr_mat(2);

    % Pearson correlation between downsampled theta amplitude and DF/F for
    % each ROI
    tmp_corr_mat                  = corrcoef(ds_theta_ampl', dff_rem{1}(roi_nr, :));
    theta_ampl_corr_ds(roi_nr, 1) = tmp_corr_mat(2);
end

% Plot results
figure, 

sgtitle('Correlation between downsampled REM theta band/amplitude with DF/F', 'FontSize', 10)

subplot(221)
histogram(theta_corr_ds, BinMethod='fd')
ylabel('Counts')
xlabel('Theta corr. coef')

subplot(222)
boxplot(theta_corr_ds);
ylabel('Theta corr. coef')

subplot(223)
histogram(theta_ampl_corr_ds, BinMethod='fd')
ylabel('Counts')
xlabel('Theta corr. coef')

subplot(224)
boxplot(theta_ampl_corr_ds);
ylabel('Theta corr. coef')

%% Strategy 2: Upsample DF/F traces to match ephys. 
upsample_factor = ephys_srate/imaging_srate;

pnts = length(dff_rem{1}(1,:));

new_pnts = pnts*upsample_factor;

new_time = (0:new_pnts-1)/(upsample_factor*imaging_srate);

new_time(new_time>time_imaging{1}(end)) = [];

new_srate_actual = 1/mean(diff(new_time));

% Preallocate
[theta_corr_us, theta_ampl_corr_us] = deal( zeros( size(dff_rem{1},1), 1));

% Loop over ROIs and compute theta correlation
for roi_nr = 1:size(dff_rem{1},1)
    
    % Create interpolation object
    F = griddedInterpolant(time_imaging{1}, dff_rem{1}(roi_nr,:), 'spline');
    
    % Interpolate/upsample
    upsampled_dff = F(new_time);

    % Find difference between upsampled DF/F and ephys. Usually few samples.
    % Subtract them from ephys so that signal lengths match. 
    signal_diff = length(rem_theta_lfp{1}) - length(upsampled_dff);

    % Correlate upsampled DF/F and theta
    tmp_corr_mat = corrcoef( upsampled_dff, rem_theta_lfp{1}(1:end-signal_diff));
    theta_corr_us(roi_nr, 1) = tmp_corr_mat(2);

    % Correlate upsampled DF/F and theta amplitude
    tmp_corr_mat = corrcoef( upsampled_dff, rem_theta_ampl{1}(1:end-signal_diff));
    theta_ampl_corr_us(roi_nr, 1) = tmp_corr_mat(2);

end

% Plot results
figure, 

sgtitle('Correlation between upsampled DF/F and REM theta band/amplitude', 'FontSize', 10)

subplot(221)
histogram(theta_corr_us, BinMethod='fd')
ylabel('Counts')
xlabel('Theta corr. coef')

subplot(222)
boxplot(theta_corr_us);
ylabel('Theta corr. coef')

subplot(223)
histogram(theta_ampl_corr_us, BinMethod='fd')
ylabel('Counts')
xlabel('Theta corr. coef')

subplot(224)
boxplot(theta_ampl_corr_us);
ylabel('Theta corr. coef')

% Sort correlations (theta band & theta band amplitude envelope)
[sorted_corr, sorted_corr_idx]           = sort(theta_corr_us, 'descend');
[sorted_corr_ampl, sorted_corr_idx_ampl] = sort(theta_ampl_corr_us, 'descend');