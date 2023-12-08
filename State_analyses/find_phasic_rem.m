function output = find_phasic_rem(sData)

%{
Find phasic REM events in REM sleep based on procedure in de Almeida-Filho
et al. (2021).
%}


%% Step 1: Get raw and theta filtered LFP

lfp       = sData.ephysdata.lfp;
lfp_theta = sData.ephysdata.thetaband;
all_rem_start_stop = rem_sleep(sData);

lfp_delta         = sData.ephysdata.deltaband;
theta_ampl        = abs(hilbert(lfp_theta));
delta_ampl        = abs(hilbert(lfp_delta));
theta_ampl_smooth = smoothdata(theta_ampl, 'gaussian', 5000);
delta_ampl_smooth = smoothdata(delta_ampl, 'gaussian', 5000);
theta_Delta_R     = theta_ampl_smooth./delta_ampl_smooth;

for n_rem_ep = 1:size(all_rem_start_stop(:,1),1)

    rem_start_stop = all_rem_start_stop(n_rem_ep,1):all_rem_start_stop(n_rem_ep,2);

    thetaR_rem = theta_Delta_R(rem_start_stop);
    thetaA_rem = theta_ampl_smooth(rem_start_stop);
    
    lfp_theta_rem = lfp_theta(rem_start_stop);
    lfp_rem       = lfp(rem_start_stop);
    
    time_vec_rem = (0:length(lfp_theta_rem)-1)./2500;
    time_vec     = (0:length(lfp_theta)-1)./2500;
    
    %% Step 2: Find theta peaks via positive-to-negative zero crossings in the
    % derivative of filtered data
    
    theta_diff  = diff(lfp_theta_rem);
    theta_peaks = find(theta_diff(1:end-1) >= 0 & theta_diff(2:end) < 0);    % Returns Down Zero-Crossing Indices
    
    %% Step 3: Smooth inter-peak interval
    kernelSize        = 11;
    rectangularKernel = ones(1, kernelSize);
    rectangularKernel = rectangularKernel / sum(rectangularKernel); 
    
    interPeakIntervals = diff(theta_peaks);
    smoothedInterPeakIntervals = zeros(size(interPeakIntervals));
    
    % Apply smoothing to each inter-peak interval using the rectangular kernel
    for i = 1:length(interPeakIntervals)
        startIndex = max(1, i - floor(kernelSize/2));
        endIndex   = min(length(interPeakIntervals), i + floor(kernelSize/2));
        
        % Apply the rectangular kernel to the current inter-peak interval
        smoothedInterPeakIntervals(i) = dot(rectangularKernel(1:endIndex-startIndex+1),  interPeakIntervals(startIndex:endIndex));
    
    end
    
    % manual smoothing (convolution):
    
    nSign = length(interPeakIntervals);
    nKern = length(rectangularKernel);
    nConv = nSign + nKern - 1;
    
    half_kern = floor(nKern/2);
    % kflip     = rectangularKernel(end:-1:1);
    % data4conv = [zeros(1, half_kern) interPeakIntervals' zeros(1, half_kern)];
    % conv_res = zeros(1, nConv);
    % 
    % for ti = half_kern+1:nConv-half_kern
    % 
    %     tmpdata = data4conv(ti-half_kern:ti+half_kern);
    % 
    %     conv_res(ti) = sum( tmpdata .* kflip);
    % end
    % 
    % dataX   = fft(interPeakIntervals, nConv);  
    % kernelX = fft(kflip, nConv);
    % 
    % convres = ifft( dataX .* kernelX);
    % convres = convres(nKern+1:end-nKern);
    
    
    smoothedInterPeakIntervalsM = conv(rectangularKernel, interPeakIntervals, 'full');
    smoothedInterPeakIntervalsM = smoothedInterPeakIntervalsM(half_kern+1:end-half_kern);
    %% Plot inter-peak intervals and histogram
    % figure;
    % h(2) = subplot(2,1,1);
    % plot(interPeakIntervals, '-o');
    % title('Original Inter-Peak Intervals');
    % xlabel('Peak Index');
    % ylabel('Inter-Peak Interval');
    % 
    % h(1) = subplot(2,1,2); hold on
    % plot(smoothedInterPeakIntervals, '-o');
    % yline(prctile(smoothedInterPeakIntervals, 10), 'r--', 'LineWidth',1)
    % 
    % % plot(n_theta_peaks)
    % title('Smoothed Inter-Peak Intervals');
    % xlabel('Peak Index');
    % ylabel('Smoothed Inter-Peak Interval');
    % linkaxes(h, 'x')
    % 
    % figure,
    % histogram(smoothedInterPeakIntervals)
    % xline(prctile(smoothedInterPeakIntervals, 10), 'r--', 'LineWidth',1)
    % xlabel('Smoothed inter-peak intervals')
    
    %% Find inter-peak intervals below 10th percentile of distribution
    data_below10_idx = smoothedInterPeakIntervals < prctile(smoothedInterPeakIntervals, 10);
    
    
    [putative_phREM_start, putative_phREM_stop] = findTransitions(data_below10_idx);
    
    n_theta_peaks = 1:length(data_below10_idx);
    
    % Preallocate
    [tmp_peak_idx, tmp_theta_peak_ephys_locs, ]                            = deal ( cell(1, length(putative_phREM_start)));
    [tmp_putative_phREM_durations, tmp_min, tmp_putative_phREM_mean_theta] = deal( zeros(1, length(putative_phREM_start)));
    
    for i = 1:length(putative_phREM_start)
    
        tmp_peak_idx{1,i} = n_theta_peaks(putative_phREM_start(i):putative_phREM_stop(i) );
    
        tmp_theta_peak_ephys_locs{1,i}    = theta_peaks(tmp_peak_idx{1,i});
        tmp_putative_phREM_durations(i) = tmp_theta_peak_ephys_locs{1,i}(end)-tmp_theta_peak_ephys_locs{1,i}(1);
    
        tmp_min(i) = min( smoothedInterPeakIntervals(tmp_peak_idx{1,i}));
    
        tmp_putative_phREM_mean_theta(i) = mean( abs( hilbert( lfp_theta_rem( tmp_theta_peak_ephys_locs{1,i}(1):tmp_theta_peak_ephys_locs{1,i}(end) ))));
    
    end
    
    %% Apply heuristics 
    
    % Heuristic # 1: duration longer than 900ms
    heuristic1_idx = tmp_putative_phREM_durations./2500 > 0.9;
    
    % Heuristic # 2: minimum inter-peak interval within candidate phREM shorter
    % than the 5th percentile of smoothed inter-peak intervals
    
    fifth_perc_threshold = prctile( smoothedInterPeakIntervals, 5);
    heuristic2_idx       = tmp_min < fifth_perc_threshold;
    
    % Heuristic # 3: mean amplitude within theta band of candidate phREM higher
    % than the mean theta amplitude of all REM sleep
    mean_rem_theta_ampl = mean( abs(hilbert(lfp_theta_rem)));
    heuristic3_idx      = tmp_putative_phREM_mean_theta > mean_rem_theta_ampl;
    
    % Phasic REM eevents
    tmp_idx      = vertcat( heuristic1_idx, heuristic2_idx, heuristic3_idx);
    phREM_idx    = sum(tmp_idx,1) == 3;
    phREM_events = tmp_theta_peak_ephys_locs(phREM_idx);
    
    
    [lfp_theta_rem_phasic, lfp_rem_phasic] = deal( NaN(1,length(lfp_theta_rem)));
    
    if isempty(phREM_events)
        msgbox('No phasic REM events found!')
        continue
    
    end
    
    for j = 1:length(phREM_events)
        
        tmp_snippet = phREM_events{1,j}(1):phREM_events{1,j}(end);
    
        lfp_theta_rem_phasic( tmp_snippet)   = lfp_theta_rem(tmp_snippet );
        lfp_rem_phasic(  tmp_snippet) = lfp_rem( tmp_snippet );
    
        phREM_lengths_seconds(j) = length(tmp_snippet(1):tmp_snippet(end))./2500;
    end
    
    total_rem_length_seconds = length( rem_start_stop )./2500;
    
    phREM_percentage = (sum( phREM_lengths_seconds(:)) / total_rem_length_seconds) * 100;
    
    %% Plot phasic REM epochs over entire REM episode
    figure, 
    h(1) = subplot(211); 
    plot(time_vec_rem, lfp_rem, 'k'), hold on
    plot(time_vec_rem, lfp_rem_phasic, 'r')
    title(sData.sessionInfo.sessionID, 'interpreter', 'none')
    ylabel('Raw LFP')
    
    h(2) = subplot(212);
    plot(time_vec_rem, lfp_theta_rem, 'k'), hold on
    plot(time_vec_rem, lfp_theta_rem_phasic, 'r')
    plot(time_vec_rem, thetaR_rem*0.1-1 )
    plot(time_vec_rem, thetaA_rem*6-1)
    xlabel('Time from REM start (s)', FontSize=16)
    ylabel('Theta (5-12Hz) filtered LFP')
    legend('Filtered LFP','Phasic REM events', 'theta/delta Ratio', 'Theta ampl.')
    title(['Phasic REM = ', num2str(phREM_percentage), ' % of total REM']);
    linkaxes(h, 'x')
    
    %% Optionally, if you want to obtain the smoothed signal with peaks included
    % smoothedSignal = zeros(size(lfp_theta_rem));
    % smoothedSignal(theta_peaks(data_below10_idx)) = .5;
    % 
    % figure, 
    % % h(1) = subplot(211);
    % % plot(time_vec_rem, lfp_rem), hold on
    % plot(time_vec_rem, lfp_theta_rem), hold on
    % 
    % plot(time_vec_rem, smoothedSignal)
    % 
    % h(2) = subplot(212);
    % plot(time_vec_rem, lfp_theta_rem), hold on
    % plot(time_vec_rem, smoothedSignal)
    
    % linkaxes(h, 'x')
    % test = zeros(1, length(theta_diff));
    % test(theta_peaks) = .1;
end


