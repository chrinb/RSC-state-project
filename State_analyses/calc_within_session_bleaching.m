function [mean_dff_beginning, mean_dff_end, mean_dff, per_cell_dff_beginning, per_cell_dff_end] = calc_within_session_bleaching(sData)

%{
Investigate decrease in fluorescence across the imaging session
%}

dff = sData.imdata.roiSignals(2).newdff;

mean_dff = mean(dff, 1, 'omitnan');
n_frames = 200;

imaging_sampling_rate = find_imaging_framerate(sData);

%% Compare the mean DF/F in the first 200 frames with last 200 frames 

mean_dff_beginning = mean( mean_dff( 1:n_frames), 'omitnan');
mean_dff_end       = mean( mean_dff(end-n_frames:end), 'omitnan');

per_cell_dff_beginning = mean(dff(:, 1:n_frames), 2, 'omitnan');
per_cell_dff_end       = mean(dff(:, end-n_frames:end),2,'omitnan');

%% Plot mean DF/F
if isfield(sData, 'episodes')
    REM_episodes = rem_sleep(sData);
    NREM_episodes = nrem_sleep(sData);

    NREM_start_end = round( (NREM_episodes./2500)*imaging_sampling_rate);
    REM_start_end = round( (REM_episodes./2500)*imaging_sampling_rate);
end


figure, 
sgtitle(sData.sessionInfo.sessionID, 'interpreter', 'none')
subplot(211), hold on
plot( mean( dff, 'omitnan'));

y_max = max(mean(dff, 'omitnan'));
y_min = min(mean(dff, 'omitnan'));

 if isfield(sData, 'episodes')
        % NREM bouts
        for i = 1:length(NREM_start_end)
            x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
            y = [y_min y_max y_max y_min];
            patch(x, y, 'blue', 'edgecolor', 'none', 'FaceAlpha', .2);
        end
        % REM bouts
        if size(REM_start_end(:,1),1) == 1
            a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)];
            b = [y_min y_max y_max y_min];
            patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
        else
            for i = 1:length(REM_start_end)
                a = [REM_start_end(i,1) REM_start_end(i,1) REM_start_end(i,2) REM_start_end(i,2)];
                b = [y_min y_max y_max y_min];
                patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
            end
        end
 end
xlabel('Frames', FontSize=14)
ylabel('Mean DF/F', FontSize=14)
%% Plot binned mean DF/F
tmp_data_imag_mean = mean(dff,'omitnan');

bin_win_sec                                  = 6;
bin_win_frames                               = bin_win_sec*round(imaging_sampling_rate);
pnts_to_remove_imag                          = mod( size(mean(dff),2), bin_win_frames);
data_trim_imag                               = tmp_data_imag_mean(:, 1:end-pnts_to_remove_imag);
data_trim_imag                               = tmp_data_imag_mean(:, 1:end-pnts_to_remove_imag);
reshaped_data                                = reshape(data_trim_imag,  bin_win_frames, []);
binned_data                                  = mean( reshaped_data,1);

subplot(212)
plot( binned_data);
xlabel(['Binned frames (bin size = ' num2str(bin_win_sec), ' sec)'], FontSize=14)
ylabel('Binned Mean DF/F', FontSize=14)