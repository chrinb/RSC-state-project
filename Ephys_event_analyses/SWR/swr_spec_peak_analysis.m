figure(1), clf 

ecog_pre = mean_data_ecog_pre;
ecog_post = mean_data_ecog_post;
lfp_pre   = mean_data_lfp_pre;
lfp_post = mean_data_lfp_post;
% time = time;

findpeaks(lfp_pre, time,  'MinPeakHeight', max(lfp_pre)*0.99);
[lfp_run_peak, lfp_pre_loc] = findpeaks(lfp_pre, time, 'MinPeakHeight', max(lfp_pre)*0.99);
hold on

findpeaks(ecog_pre,time, 'MinPeakHeight', max(ecog_pre)*0.99);
[ecog_run_peak, ecog_pre_loc] = findpeaks(ecog_pre, time, 'MinPeakHeight', max(ecog_pre)*0.99);

findpeaks(lfp_post, time,  'MinPeakHeight', max(lfp_post)*0.99);
[lfp_sleep_peak, lfp_post_loc] = findpeaks(lfp_post, time, 'MinPeakHeight', max(lfp_post)*0.99);

findpeaks(ecog_post, time,'MinPeakHeight', max(ecog_post)*0.99);
[ecog_sleep_peak, ecog_post_loc] = findpeaks(ecog_post, time, 'MinPeakHeight', max(ecog_post)*0.99);


titleString = sprintf(['LFP pre = ', num2str(lfp_pre_loc), ', ECoG pre = ', num2str(ecog_pre_loc), ...
    ', lfp post = \n',  num2str(lfp_post_loc), ', ecog post = ', num2str(ecog_post_loc)]);
title(titleString)

% ecog_sleep_loc-lfp_sleep_loc
ecog_pre_loc-lfp_pre_loc
ecog_post_loc-lfp_post_loc
xlabel('Time from SWR peak (s)')
set(gca, 'FontSize',10)
legend({'LFP pre','', 'ECoG pre','', 'lfp post','', 'ecog post'})

%% Unit scale

lfp_pre_scaled  = (lfp_pre - min(lfp_pre) ) / ( max(lfp_pre) - min(lfp_pre));
ecog_pre_scaled = (ecog_pre - min(ecog_pre) ) / ( max(ecog_pre) - min(ecog_pre));

lfp_post_scaled  = (lfp_post - min(lfp_post) ) / ( max(lfp_post) - min(lfp_post));
ecog_post_scaled  = (ecog_post - min(ecog_post) ) / ( max(ecog_post) - min(ecog_post));

ecog_scaled = (ecogtm - min(ecogtm) ) / ( max(ecogtm) - min(ecogtm));
lfp_scaled  = (lfptm - min(lfptm) ) / ( max(lfptm) - min(lfptm));

figure(2), clf
hold on
plot(time, lfp_pre_scaled, 'r')
xline(lfp_pre_loc, 'r')
plot(time, ecog_pre_scaled, 'b')
xline(ecog_pre_loc, 'b')
plot(time, lfp_post_scaled, 'm')
xline(lfp_post_loc, 'm')
plot(time, ecog_post_scaled, 'g')
xline(ecog_post_loc, 'g')
set(gca, 'xlim', [-.05 .05], 'YLim',[0 1.1], 'FontSize', 14)
legend({'LFP pre','', 'ECoG pre','', 'lfp post','', 'ecog post',''})
xlabel('Time from SWR peak (s)')
ylabel('150-300Hz mean power (normalized)')
title('Unity-normed data ')

%% Visualize some mean ripple band signals
figure(10)
for i = 1:size(ecogt,3)
    plot(time, ecogt(1,:,i))
    hold on
    plot(time, lfpt(1,:,i))
    title( num2str(i))
    pause(.7)
    clf
end


%% Find peaks

% Average over frequencies 150 Hz and above
lfp = squeeze( mean( all_specs_lfp(6:end, :,:), 1));
ecog = squeeze( mean( all_specs_ecog(6:end, :,:), 1));

% Time indices for -0.05 to 0.05 s after SWR peak
[~, min_idx] = min( abs(time - -.05));
[~, min_idx2] = min( abs(time - .05));

% Preallocate
[all_lfp_peaks, all_ecog_peaks] = deal( zeros( 1, size(ecogt,3)));

figure(11), clf

% Loop over nr of peri-SWR events
for i = 1:size(ecogt,3)
    
    % Trim down LFP/ECoG snippets to -0.1 to +0.1 s
    lfp_trim = lfp( min_idx:min_idx2, i);
    ecog_trim = ecog( min_idx:min_idx2, i);
    
    % Search for peaks in interval and choose largest
    [tmp_lfp_peaks, tmp_lfp_loc] = findpeaks(lfp_trim, time(min_idx:min_idx2) );
    [~, max_id1]     = max(tmp_lfp_peaks);
    all_lfp_peaks(i) = tmp_lfp_loc(max_id1);

    % Repeat for ECoG signal
    [tmp_ecog_peaks, tmp_ecog_loc] = findpeaks(ecog_trim, time(min_idx:min_idx2) );
    [~, max_id2]      = max(tmp_ecog_peaks);
    all_ecog_peaks(i) = tmp_ecog_loc(max_id2);


    % Optional plot for visualization 
    % plot( time(min_idx:min_idx2), lfp_trim, 'k')
    % xline(tmp_lfp_loc(max_id1), 'k', 'LineWidth',1)
    % hold on
    % 
    % plot( time(min_idx:min_idx2), ecog_trim, 'm')
    % xline(tmp_ecog_loc(max_id2), 'm', 'LineWidth',1)
    % hold on
    % title( num2str(i))

    % pause(.5)
    % clf
end

%% Plot histograms
figure(12), clf

h = histogram(all_lfp_peaks, 'BinMethod','fd', 'FaceColor','r');
xline( mean(all_lfp_peaks), 'r', 'LineWidth',1)
hold on

k = histogram(all_ecog_peaks, 'BinMethod','fd', 'FaceColor','b');
k.BinEdges = h.BinEdges;

xline( mean(all_ecog_peaks), 'b','LineWidth',1)


[p, h] = ttest(all_lfp_peaks, all_ecog_peaks);

titlestring = sprintf( ['Mean LFP peak = ', num2str(mean(all_lfp_peaks) ), ', mean ECoG peaks = ', num2str( mean(all_ecog_peaks) ),'\n p = ', num2str(h), ' (paired t-test)' ] );
legend({'LFP peaks','' ,'ECoG peak', ''})
xlabel('Time from SWR peak')
ylabel('Counts')
set(gca, 'FontSize', 15)
title(titlestring, FontSize=12)
