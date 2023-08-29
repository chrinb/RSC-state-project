function [mean_pow_mat, mean_dff_mat] = binned_theta_dff_analysis(sData, params)

% Written by Christoffer Berge | Vervaeke lab

% Extract REM theta power signal and bin signal in 4 second bins and compute 
% mean power. Extract the corresponding DF/F data during rem, average and
% bin signal. 

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

% Get 2P frame rate
imaging_sampling_rate = find_imaging_framerate(sData);

if unique(sData.imdata.plane_indices) > 1

    % Divide by 4 to get FPS per plane
    imaging_sampling_rate = imaging_sampling_rate/4;
end

%% Bin REM theta power and DF/F and compute bin means

% Create bin windows
bin_length_sec = 2;

bin_win_imag  = round(imaging_sampling_rate*bin_length_sec);
bin_win_ephys = 2500*bin_length_sec;

params.state = 'REM';
rem_data       = freq_band_analysis(sData, params);

mean_pow_cell = cell( size(rem_data.state_lfp_filt_ampl, 1), 1);
mean_dff_cell = cell( size(rem_data.state_lfp_filt_ampl, 1), 1);

% Loop over episodes
for ep_nr = 1:size(rem_data.state_lfp_filt_ampl, 1)
    
    % Get data
    tmp_data_ephys = rem_data.state_lfp_filt_ampl{ep_nr, 1};
    
    % Continue with analysis only if current REM data is NOT empty (if it
    % is, the episode was skipped because it didn't fit the criteria).
    if ~isempty(tmp_data_ephys)

        % Start with ephys data

        % Find remainder after splitting data into bins
        pnts_to_remove_ephys = mod( size(tmp_data_ephys,1), bin_win_ephys);
        
        % Remove remainder from data
        data_trim_ephys = tmp_data_ephys(1:end-pnts_to_remove_ephys);
        
        % Split data into bins
        reshaped_data_ephys = reshape(data_trim_ephys,  bin_win_ephys, []  );
    
        % Compute mean power per bin
        mean_pow = mean( reshaped_data_ephys,1);
    
        % Add remainder
%         mean_pow_cell{ep_nr, 1} = [mean_pow, mean(tmp_data_ephys(end-pnts_to_remove_ephys:end))];
        mean_pow_cell{ep_nr, 1} = mean_pow;

        % Now for imaging data

        % Get start/end of REM in imaging time
        rem_snippet_imaging = frames(rem_data.state_times(ep_nr, 1)):frames(rem_data.state_times(ep_nr, 2));
        
        % Extract imaging data during REM episode and compute mean over ROIs
        tmp_data_imag = dff(:, rem_snippet_imaging);
        tmp_data_imag = mean(tmp_data_imag, 1, 'omitnan');
        
        % Find remainder after splitting data into bins
        pnts_to_remove_imag = mod( size(tmp_data_imag,2), bin_win_imag);
        
        % Remove remainder from data
        data_trim_imag = tmp_data_imag(:, 1:end-pnts_to_remove_imag);
        
        % Split data into bins
        reshaped_data_imag = reshape(data_trim_imag,  bin_win_imag, []);
    
        % Compute mean power per bin
        mean_dff = mean( reshaped_data_imag,1);
    
        % Add remainder
%         mean_dff_cell{ep_nr, 1} = [mean_dff, mean(tmp_data_imag(end-pnts_to_remove_imag:end))];
        mean_dff_cell{ep_nr, 1} = mean_dff;

        % Because of different sampling rates, the nr of bins may differ despite
        % similar bin window length. Therefore find the dataset with fewest bins,
        % and throw out additional bin data in the larger dataset. 
        min_bin_nr = min( [size(reshaped_data_imag, 2), size(reshaped_data_ephys, 2)]);

        mean_pow_cell = mean_pow_cell{ep_nr, 1}(1:min_bin_nr);
        mean_dff_cell = mean_dff_cell{ep_nr, 1}(1:min_bin_nr);

    end

end

% Convert from cell array to matrix
mean_pow_mat = cell2mat( mean_pow_cell);
mean_dff_mat = cell2mat( mean_dff_cell);


%% Calculate linear regression
% x = mean_pow_mat{1}';
% y = mean_dff_mat{1}';
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