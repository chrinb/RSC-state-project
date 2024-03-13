function [power_lfp, freq_lfp, power_ecog, freq_ecog] = plot_qw_power(sData, params)


if strcmp(params.qw_or_run, 'qw')
    epochs  = ~(sData.analysis.movement_threshold);
elseif strcmp(params.qw_or_run, 'run')
    epochs = sData.analysis.paw_threshold;
end

lfp  = sData.ephysdata.lfp;
ecog = sData.ephysdata2.lfp;

% Z-score to remove DC offset
lfp  = zscore(lfp);
ecog = zscore(ecog);

% Convert to ephys time
[paw_start, paw_stop] = findTransitions(epochs);

frames = sData.daqdata.frame_onset_reference_frame;

ephys_srate = 2500;

for i = 1:length(paw_start)

    % Find closest matching ephys
    [~, min_idx_start(i)] = min( abs( paw_start(i) - frames  ));
    [~, min_idx_stop(i)] = min( abs( paw_stop(i) - frames  ));

    lfp_snip{i}  = lfp(min_idx_start(i):min_idx_stop(i)); 
    ecog_snip{i} = ecog(min_idx_start(i):min_idx_stop(i));

    % Skip too short segments to perform pwelch with default window and overlap 
    if length(lfp_snip{i}) > 8
        [power_lfp{i}, freq_lfp{i}]   = pwelch( lfp_snip{i}, [],[],[], ephys_srate);
        [power_ecog{i}, freq_ecog{i}] = pwelch( ecog_snip{i}, [],[],[], ephys_srate);
    else
        [power_lfp{i}, freq_lfp{i}]   = deal([]);
        [power_ecog{i}, freq_ecog{i}] = deal([]);
    end
end

% 
% all_data = horzcat( horzcat( power_lfp{1, :}));
% all_f    = horzcat( horzcat( freq_lfp{1, :}));
% 
% % Find index and length of longest episode
% % [max_length, max_idx] = max( cellfun(@(x) length(x), all_data), [], 'all' ); 
% [max_length, max_idx] = max( cellfun(@(x) length(x), power_lfp), [], 'all' ); 
% 
% % In the corresponding max frequencies vector, find index closest to
% % 30Hz 
% max_freq_vec = freq_lfp{max_idx};
% [test, testI] = min( abs(max_freq_vec-30));
% 
% max_vec = power_lfp{max_idx};
% max_f   = freq_lfp{max_idx};
% % % Loop over states
% % for state_nr = 1:4
% %     tmp_state_data  = all_data(state_nr,:);
% %     tmp_state_freqs = all_f(state_nr,:);
% 
%     % Loop over episodes and interpolate
%     % interpolated_episodes = {data_dim(1) , data_dim(2)};
%     for i = 1:numel(power_lfp)
% 
%         tmp_power = power_lfp{i};
%         tmp_freqs = freq_lfp{i};
%         [~, tmp_idx] = min( abs(tmp_freqs-30));
% 
%         if ~isempty(tmp_power) && length(tmp_power) < max_length
%             x  = 1:tmp_idx;
%             v  = tmp_power(1:tmp_idx);
%             xq = (linspace(1, length(tmp_power(1:tmp_idx)), testI));
%             interpolated_episodes{i} = interp1(x, v, xq, 'spline' );
%         end
%     end
% % end
% 
% 
% data_cat  = vertcat( interpolated_episodes{:});
% data_mean = mean(data_cat);
% data_sem  = std(data_cat, 0,1) ./ sqrt( size(data_cat,1));
% 
% plot_colors = {'r', 'b'};

% figure, 
% % plot(data_mean)
% 
% % plot( max_f(1:testI)', mean(data_cat{1, state_nr})  )
% % shadedErrorBar( max_f(1:testI), data_mean, data_sem, 'lineProps', {'markerfacecolor', plot_colors{state_nr} } )
% shadedErrorBar( max_f(1:testI), data_mean, data_sem )
% 
% set(gca, 'xlim',[0 20])
% xlabel('Frequency (Hz)')
% ylabel('Power \muV^2')
% ax = gca;
% ax.FontSize = 16;
% axis square


%% Plot epochs in face camera vector and in ephys as sanity check
% srate = find_imaging_framerate(sData);
% cam_fps = 30.94;
% 
% time_cam     = linspace(1, length(paw_epochs), length(paw_epochs))/cam_fps;
% % time_vec = (0:length(paw_epochs)-1)/cam_fps;
% 
% 
% % time_imaging = linspace(1, numel(sData.daqdata.frameIndex), numel(sData.daqdata.frameIndex) )/srate;
% time_ephys   = linspace(1, length(lfp), length(lfp) )/2500;
% 
% figure(1), clf,
% hold on
% h(1) = subplot(211);
% plot(time_cam, paw_epochs)
% 
% [start_vec, stop_vec] = deal( zeros(1, size(time_ephys,2)));
% start_vec(min_idx_start) = 1;
% stop_vec(min_idx_stop) = 1;
% title('Face camera')
% set(gca, 'ylim', [0 2])
% 
% h(2) = subplot(212); hold on
% plot( time_ephys, start_vec, 'r')
% plot( time_ephys, stop_vec, 'b')
% title('Ephys')
% xlabel('Time (s)')
% set(gca, 'ylim', [0 2])
% legend('Start', 'Stop')
% linkaxes(h, 'x')
