function plot_peri_swr_spectrogram(spectrogram_data, params, time, pre_idx, post_idx, across_day_idx)   

%{ 
Plot (1) the average ripple band power from the spectrogram outout and its 
associated SEM; and (2) spectrograms of pre- vs post peri SWR LFP/ECoG (entire
recording traces are first z-scored) in log 10 scale. Top row is "standard"
signal, bottom have z-scored each frequency across time)
%}

 % Split data into pre and post
 pre_spec     = vertcat( spectrogram_data{pre_idx});
 post_spec    = vertcat( spectrogram_data{post_idx, 1});

 frequencies  = spectrogram_data{1, 1}{1, 2};
 n_freqs      = numel( spectrogram_data{1, 1}{1, 2} );
 n_timepoints = numel( spectrogram_data{1, 1}{1, 1} );


% Concatenate SWR trials
 [pre_spec_cat, post_spec_cat] = deal( zeros(n_freqs, n_timepoints, 0));
 for i = 1:size(pre_spec, 1)
     % Concatenate along the third dimension
     pre_spec_cat  = cat(3, pre_spec_cat, pre_spec{i, 4});
     post_spec_cat = cat(3, post_spec_cat, post_spec{i, 4});
 end

if strcmp(params.signal_type, 'ECoG' )
    fig_nr1 = 12;
    fig_nr2 = 13;
else
    fig_nr1 = 14;
    fig_nr2 = 15;
end


% Average pre and post
pre_spec_mean  = mean(pre_spec_cat, 3);
post_spec_mean = mean(post_spec_cat, 3); 


lower_ripple_freq = 150;
[~, lower_ripple_freq_idx] = min( abs( spectrogram_data{1, 1}{1, 2} - lower_ripple_freq)) ;

all_specs      = cat(3, pre_spec_cat, post_spec_cat);
% extract PSD from 150 - 300 Hz 
all_specs_mean = mean(all_specs(lower_ripple_freq_idx:end,: ,:) ,3); 
all_specs_sem   = std(all_specs(lower_ripple_freq_idx:end,: ,:),0, 3) ./ sqrt( size(all_specs(lower_ripple_freq_idx:end,: ,:),3));
            
pre_spec_mean_rip      = mean( pre_spec_mean(lower_ripple_freq_idx:end,: ,:) ,3); 
pre_spec_mean_rip_sem  = std( pre_spec_mean(lower_ripple_freq_idx:end,: ,:), 0, 3) ./ sqrt( size( pre_spec_mean(lower_ripple_freq_idx:end,: ,:),3) );
            
post_spec_mean_rip     = mean( post_spec_mean(lower_ripple_freq_idx:end,: ,:) ,3); 
post_spec_mean_rip_sem = std( post_spec_mean(lower_ripple_freq_idx:end,: ,:), 0, 3) ./ sqrt( size( post_spec_mean(lower_ripple_freq_idx:end,: ,:),3) );

figure(fig_nr1), clf
hold on
% shadedErrorBar(time, mean(pre_spec_mean_rip, 1), mean(pre_spec_mean_rip,1), 'lineprops', 'r'), 
% shadedErrorBar(time, mean(post_spec_mean_rip, 1), mean(post_spec_mean_rip,1), 'lineprops', 'b'), 
plot(time, mean(pre_spec_mean_rip), 'r')
plot(time, mean(post_spec_mean_rip), 'b')
xlabel('Time from SWR peak (s)')
set(gca, "FontSize", 12)
legend({'Pre', 'Post'})


pre_spec_mean_zscore  = zscore(pre_spec_mean, 0, 2); 
post_spec_mean_zscore = zscore(post_spec_mean, 0, 2); 

% Plot spectrogram
% time = pre_spec{1, 1};  
F = frequencies;



figure(fig_nr2), clf
% sgtitle(sessionID),
subplot(2,2,1)
contourf(time,F,10*log10(abs(pre_spec_mean)),200,'edgecolor','none');
xline(0, 'r--', 'LineWidth',2)
colormap(viridis)
colorbar
ylabel('Frequency (Hz)', FontSize=12)
title('Pre-RUN', FontSize=12)

subplot(2,2,2)
contourf(time,F,10*log10(abs(post_spec_mean)),200,'edgecolor','none');
colormap(viridis)
xline(0, 'r--', 'LineWidth',2)

colorbar
title('Post-RUN', FontSize=12)

subplot(2,2,3)
contourf(time,F,10*log10(abs(pre_spec_mean_zscore)),200,'edgecolor','none');
xline(0,'r--', 'LineWidth',2)
colormap(viridis)
caxis([-10 10])
colorbar
ylabel('z-scored frequency (Hz)', FontSize=12)
xlabel('Time from SWR peak (s)', FontSize=12)

subplot(2,2,4)
contourf(time,F,10*log10(abs(post_spec_mean_zscore)),200,'edgecolor','none');
xline(0, 'r--', 'LineWidth',2)
colormap(viridis)
caxis([-10 10])
colorbar
xlabel('Time from SWR peak (s)', FontSize=12)

%% Plot across day changes

% Sort data across days

pre_spec_data = pre_spec(:, 4);
post_spec_data = post_spec(:, 4);

% 1st and 2nd column are pre and post data from mouse 1 etc.
% pre_across_mice = reshape(pre_spec_data, 7, []);
% post_across_mice = reshape(post_spec_data, 7, []);

% pre_rf = pre_across_mice(1:4, :);
% pre_gol= pre_across_mice(5:7, :);

% post_rf  = post_across_mice(1:4, :);
% post_gol = post_across_mice(5:7, :);

% Loop over Pre and post
% for kk = 1:2

across_day_idx_trim  = across_day_idx(pre_idx);
n_days               = unique(across_day_idx_trim);

[mean_peri_rip_band_pow_per_day_pre, mean_peri_rip_band_pow_per_day_post] = deal( []);
for day_n = 1:numel(n_days)
    
    tmp_day_idx = n_days(day_n);

    day_idx      = across_day_idx_trim == tmp_day_idx;
    tmp_dat_pre  = pre_spec_data(day_idx);
    tmp_dat_post = post_spec_data(day_idx);

    % tmp_dat_pre  = pre_across_mice{day_n};
    % tmp_dat_post = post_across_mice{day_n};

    % Loop through each cell in the cell array      
    [across_mice_cat_pre, across_mice_cat_post] = deal( zeros(n_freqs, n_timepoints, 0));
    for i = 1:size(tmp_dat_pre,1)      
         
        % across_mice_cat_pre = cat(3, across_mice_cat_pre, tmp_dat_pre{day_n}{1, i}  );
        % across_mice_cat_post = cat(3, across_mice_cat_post, tmp_dat_post{1, day_n}{1, i}  );
       
        across_mice_cat_pre = cat(3, across_mice_cat_pre, tmp_dat_pre{i}  );
        across_mice_cat_post = cat(3, across_mice_cat_post, tmp_dat_post{i}  );
    
        % Concatenate along the third dimension
        % pre_spec_cat  = cat(3, pre_spec_cat, pre_spec{i, 4});
        % post_spec_cat = cat(3, post_spec_cat, post_spec{i, 4});
    end
       
         % Calculate avg across SWRs
         acr_mice_mean_pre     = mean(across_mice_cat_pre, 3);
         acr_mice_mean_pre_rip = mean(acr_mice_mean_pre(lower_ripple_freq_idx:end, :) ,1); % extract PDS from 100 - 300 Hz
         % acr_mice_sem_pre   = std(acr_mice_mean_pre(lower_ripple_freq_idx:end,: ,:),0, 3)./sqrt( size(acr_mice_mean_pre(11:end,: ,:),3));

         acr_mice_mean_post     = mean(across_mice_cat_post, 3);
         acr_mice_mean_post_rip = mean(acr_mice_mean_post(lower_ripple_freq_idx:end, :) ,1); % extract PDS from 100 - 300 Hz
         % acr_mice_sem_post   = std(acr_mice_mean_post(lower_ripple_freq_idx:end,: ,:),0, 3)./sqrt( size(acr_mice_mean_post(11:end,: ,:),3));


         mean_peri_rip_band_pow_per_day_pre(day_n, :) = acr_mice_mean_pre_rip;
         % day_cat_pre{day_n, 2} = acr_mice_sem_pre;

         mean_peri_rip_band_pow_per_day_post(day_n, :) = acr_mice_mean_post_rip;
         % day_cat_post{day_n, 2} = acr_mice_sem_post;

end

%% Plot
 figure(17), clf
for day_n = 1:numel(n_days)
    subplot(1,numel(n_days),day_n)
    hold on
    plot(time, mean_peri_rip_band_pow_per_day_pre(day_n, :), 'r'), 
    plot(time, mean_peri_rip_band_pow_per_day_post(day_n, :), 'b'), 
    % set(gca, 'ylim', [2.6e-07, 3.4e-07])
    axis square
end           
legend({'Pre', 'Post'})
sgtitle('Mean peri-SWR 150-100Hz power in RSC ECoG') 

%% sanity check: compute overall means for RF and GOL
% test_rf = vertcat(day_cat_pre{1:4,1});
% test_rf2 = vertcat(day_cat_post{1:4,1});
% 
% rf_all =[ test_rf; test_rf2];
% 
% test_gol = vertcat(day_cat_pre{5:7,1});
% test_gol2 = vertcat(day_cat_post{5:7,1});
% 
% gol_all =[ test_gol; test_gol2];
