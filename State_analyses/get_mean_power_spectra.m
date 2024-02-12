function [power_calc, frequencies_calc] = get_mean_power_spectra(sData, params)

%{
Compute mean power spectra for different states in session
%}

%%  Load variables
switch params.ephys_signal
    case 'lfp'
    signal = zscore(sData.ephysdata.lfp);
    case 'ecog'
    signal  = zscore(sData.ephysdata2.lfp);
end

srate = 2500;

ephys_state_vectors = sData.behavior;

% Skip sessions were there are large outliers (likely to be noise) 
if ~( max(signal) > std(signal)*8)
    %% Compute power spectra
    states = fieldnames(sData.behavior);
    
    window = srate*4;
    
    colors = copper(4);
    
    % Loop over states
    for state_nr = 1:4
        
        % Logical vector containing epochs of current state
        temp_state = sData.behavior.(states{state_nr});
        
        % Check that sessions contain state type
        if ~sum(temp_state) == 0
    
            [state_start, state_stop] = findTransitions(temp_state);
            
            % Compute power spectrum per episode using Welch's method
            for episode_nr = 1:size(state_start,2)
    
                tmp_ep_data = signal(state_start(episode_nr):state_stop(episode_nr));
                % tmp_ep_lfp  = lfp(state_start(episode_nr):state_stop(episode_nr));
    
                % Mean-center to remove DC offset
                % tmp_ep_ecog = tmp_ep_ecog-mean(tmp_ep_ecog);
                % tmp_ep_lfp  = tmp_ep_lfp-mean(tmp_ep_lfp);
    
                [power_calc{state_nr, episode_nr}, frequencies_calc{state_nr, episode_nr}] = pwelch(tmp_ep_data, [],[],[], srate);
              
                % [power_lfp{state_nr, episode_nr}, f_lfp{state_nr, episode_nr}]  = pwelch(tmp_ep_lfp, [],[],[], srate);
    
                % plot(tmp_pxx{state_nr, episode_nr}, 'Color', colors(state_nr,:))
            end
        else
            power_calc{state_nr, 1} = [];
            frequencies_calc{state_nr, 1}  = [];
            % f_lfp{state_nr, 1}      = [];
            % f_ecog{state_nr, 1}     = [];
        end
    
    end
else
    power_calc{4, 1} = [];
    frequencies_calc{4, 1}  = [];
    msgbox(['Skipped session ', sData.sessionInfo.sessionID])
    % f_lfp{1, 1}      = [];
    % f_ecog{1, 1}     = [];
end

% 
% tmp1 = f_lfp(3,:); % x2
% tmp2 = power_lfp(3,:); % x1
% 
% empty_idx_aw = cellfun(@isempty, tmp2);
% tmp2      = tmp2(~empty_idx_aw);
% 
% [maxV, maxI] = max( cellfun(@(x)  length(x), tmp2 ));
% 
% for i = 1:length(tmp2)
% 
%         % Interpolate short state
%         x  = 1:length(tmp2{i});
%         v  = tmp2{i};
%         xq = (linspace(1, length(tmp2{i}), maxV));
%         temp_state_interp{i} = interp1(x, v, xq, 'spline' );
% end
% 
% figure, hold on
% state_means = mean( vertcat(temp_state_interp{:}));
% cellfun(@(x) plot(f_lfp{3,maxI}, x), mean(temp_state_interp))
% 
% figure, 
% plot(f_lfp{3,maxI}, state_means)

% set(gca, 'xlim',[0 30])
% xlabel('Frequency (Hz)')
% ylabel('Power')
%% Alternative 2: Power spectrum per episode
% nrem_episodes = nrem_sleep(sData);
% 
% 
% figure, hold on
% for episode_nr = 1:size(nrem_episodes,1)
%     tmp_ep = ecog(nrem_episodes(episode_nr,1):nrem_episodes(episode_nr,2));
% 
%     tmp_pxx1{episode_nr} = pwelch(tmp_ep,  [],[],[], srate);
% 
%     plot(tmp_pxx1{episode_nr})
% end
% 
% set(gca, 'xlim',[0 30])
% xlabel('Frequency (Hz)')
% ylabel('Power')
% %% Alternative 1: Concatenate all episodes of given state, find power spectrum
% % nrem_cat = ecog(ephys_state_vectors.NREM_vector);
% % 
% % [pxx, f] = pwelch(nrem_cat, srate);
% % 
% % plot(pxx,'k', 'LineWidth',1)
% % set(gca, 'xlim',[0 30])
% % xlabel('Frequency (Hz)')
% % ylabel('Power')
% % ax = get(gca);
% % ax.FontSize = 6;
