function [mean_pow_vec, mean_dff_vec, mean_dff_mat] = binned_theta_dff_analysis(sData, params)

% Written by Christoffer Berge | Vervaeke lab

%{
Extract frequency band power signal (delta/theta/sigma), bin signal and compute 
mean power per bin. Extract the corresponding average DF/F data during
episode, bin signal, and compute mean per bin. 
%}

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

if strcmp(params.zscore, 'yes')
    dff = zscore(dff, 0, 2);
end

% Vector to assign ephys time points to imaging time 
frames = sData.daqdata.frame_onset_reference_frame;

% Get 2P frame rate
imaging_sampling_rate = find_imaging_framerate(sData);

if unique(sData.imdata.plane_indices) > 1

    % Divide by 4 to get FPS per plane
    imaging_sampling_rate = imaging_sampling_rate/4;
end

%% Bin frequency band power and DF/F and compute mean per bin

% Select frequency band for mean power normalization
switch params.freq_band
    case 'theta'
        freq = [5, 9];
    case 'SO'
        freq = [0.1, 1];
    case 'delta'
        freq = [.5 4];
    case 'sigma'
        freq = [8, 18];
end

% Create bin windows
bin_length_sec = params.bin_win;

bin_win_imag  = round(imaging_sampling_rate*bin_length_sec);
bin_win_ephys = 2500*bin_length_sec;

% Get state/frequency band-specific data
state_data = freq_band_analysis(sData, params);

% Text for selecting data
txt = 'state_ephys_filt_ampl';

% Preallocate
[mean_pow_cell, mean_binned_dff_cell, mean_dff] = deal( cell( size(state_data.(txt), 1), 1));

% Loop over episodes
for ep_nr = 1:size(state_data.(txt), 1)
    
    % Get data
    tmp_data_ephys = state_data.(txt){ep_nr, 1};
    
    % Mean power of this frequency band during session
    freq_band_mean_pow = bandpower(state_data.signal, 2500, freq );

    % Check if state data is empty (indicating that episode didn't fit
    % criteria, or if DF/F is empty (cell type not present). 
    if ~isempty(tmp_data_ephys) && ~isempty(dff)

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

        % Get start/end of episode in imaging time
        state_snippet_imaging = frames(state_data.state_times(ep_nr, 1)):frames(state_data.state_times(ep_nr, 2));
        
        % Extract imaging data during episode and compute mean over ROIs
        tmp_data_imag      = dff(:, state_snippet_imaging);
        tmp_data_imag_mean = mean(tmp_data_imag, 1, 'omitnan');
        
        mean_dff{ep_nr, 1} = tmp_data_imag_mean;
        % Find remainder after splitting data into bins
        pnts_to_remove_imag = mod( size(tmp_data_imag_mean,2), bin_win_imag);
        
        % Remove remainder from data
        data_trim_imag = tmp_data_imag_mean(:, 1:end-pnts_to_remove_imag);
        
        % Split data into bins
        reshaped_data_imag = reshape(data_trim_imag,  bin_win_imag, []);
    
        % Compute mean power per bin
        mean_binned_dff = mean( reshaped_data_imag,1);
    
        % Add remainder
%         mean_dff_cell{ep_nr, 1} = [mean_dff, mean(tmp_data_imag(end-pnts_to_remove_imag:end))];
        mean_binned_dff_cell{ep_nr, 1} = mean_binned_dff;

        % Because of different sampling rates, the nr of bins may differ despite
        % "identical" bin window length. Therefore find the dataset with fewest bins,
        % and throw out additional bin data in the larger dataset. 
        min_bin_nr = min( [size(reshaped_data_imag, 2), size(reshaped_data_ephys, 2)]);

        mean_pow_cell{ep_nr, 1} = mean_pow_cell{ep_nr, 1}(1:min_bin_nr);
        mean_binned_dff_cell{ep_nr, 1} = mean_binned_dff_cell{ep_nr, 1}(1:min_bin_nr);

        % Normalize to mean power in this frequency band during session
        % (regardless of behavioral state)
        mean_pow_cell{ep_nr, 1} = mean_pow_cell{ep_nr, 1} ./ freq_band_mean_pow;

    end

end

% Concatenate mean binned power and DF/F values into vector shape
mean_pow_vec = horzcat( mean_pow_cell{:});
mean_dff_vec = horzcat( mean_binned_dff_cell{:});
mean_dff_mat =  mean_dff;
