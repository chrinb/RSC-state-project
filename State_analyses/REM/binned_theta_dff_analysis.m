function output = binned_theta_dff_analysis(sData, params)

% Written by Christoffer Berge | Vervaeke lab

% Extract REM theta signal (filtered and fitlered amplitude signal) and bin
% signal in 4 second bins and compute mean power. Bin corresponding mean
% DF/F in 4s bins and compute mean DF/F for those bins

%%  Load DF/F 
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

% Vector to assign ephys time points to imaging time 
frames = sData.daqdata.frame_onset_reference_frame;

%% Bin REM theta power and compute mean over bins

params.state = 'REM';
rem_data       = freq_band_analysis(sData, params);

% Bin window = 4s
bin_win = 2500*4;

% Loop over episodes
for ep_nr = 1:size(rem_data.state_lfp_filt_ampl,2)
    
    % Get data
    tmp_data_ephys = rem_data.state_lfp_filt_ampl{ep_nr, 1};
    
    % Find remainder after splitting data into bins
    pnts_to_remove = mod( size(tmp_data_ephys,1),bin_win);
    
    % Remove remainder from data
    data_trim = tmp_data_ephys(1:end-pnts_to_remove);
    
    % Split data into bins
    reshaped_data = reshape(data_trim,  bin_win, []  );

    % Compute mean power per bin
    mean_pow = mean( reshaped_data,1);

    % Add remainder
    mean_pow = [mean_pow, mean(tmp_data_ephys(end-pnts_to_remove:end))];

end

%% Bin DF/F and compute mean over bins

% Get 2P frame rate
imaging_sampling_rate = find_imaging_framerate(sData);

if unique(sData.imdata.plane_indices) > 1

    % Divide by 4 to get FPS per plane
    imaging_sampling_rate = imaging_sampling_rate/4;
end

bin_win_imag = round(imaging_sampling_rate*4);

% Loop over REM episodes
for rem_nr = 1:size(rem_data.state_times,1)
    
    % Get start/end of REM in imaging time
    rem_snippet_imaging = frames(rem_data.state_times(rem_nr, 1)):frames(rem_data.state_times(rem_nr, 2));
    
    % Extract imaging data during REM episode and compute mean over ROIs
    tmp_data_imag = dff(:, rem_snippet_imaging);
    tmp_data_imag = mean(tmp_data_imag, 1, 'omitnan');
    
    % Find remainder after splitting data into bins
    pnts_to_remove = mod( size(tmp_data_imag,2), bin_win_imag);
    
    % Remove remainder from data
    data_trim = tmp_data_imag(:, 1:end-pnts_to_remove);
    
    % Split data into bins
    reshaped_data = reshape(data_trim,  bin_win_imag, []);

    % Compute mean power per bin
    mean_dff = mean( reshaped_data,1);

    % Add remainder
    mean_dff = [mean_dff, mean(tmp_data_imag(end-pnts_to_remove:end))];
end

% Save results in struct
output = struct();

output.mean_binned_theta_pow = mean_pow';
output.mean_binned_dff       = mean_dff';

%% Calculate linear regression
% x = mean_pow';
% y = mean_dff';
% 
% X = [ ones(length(x),1 ), x];
% 
% b = X\y;
% 
% yCalc = X*b;
% 
% figure, 
% scatter(x, y)
% hold on
% plot(x, yCalc, '-')
% xlabel('Binned theta power', FontSize=16)
% ylabel('Binned mean DF/F', FontSize=16)