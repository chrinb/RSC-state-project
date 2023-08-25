function output = oscillation_coupling(varargin)

% Function that looks for sequential coupling of SWRs, delta waves/SOs, and
% sleep spindles. Analysis is based on the method outlined in Maingret et
% al. (2016) Nat Neuro. 


sData                        = varargin{1,1};
ecog                       = sData.ephysdata2.lfp;
ecog_delta                  = sData.ephysdata2.deltaband;
hpc_lfp                      = sData.ephysdata.lfp;
time                         = (0:length(ecog)-1)/2500;
srate                        = 2500;
% [nrem_SOs, nrem_delta_waves] = slow_wave_an(sData);
[nrem_SO, nrem_delta_waves] = mark_slow_wave(sData);

ECoG_spindle                 = sData.ephysdata2.NREMspindleStartEnd;
[swr_start_stop,~]           = mark_ripple_onset_offset(sData);

% if length(varargin) > 1
%     % Specify which sequence to plot
%     prompt = sprintf('Sequence to plot: delta-spindle(1) swr-delta (2) swr-delta-spindle(3)');
%     seq_to_plot = input(prompt);
% end

% Concatenate SO and delta waves for analysis
slow_wave = vertcat(nrem_SO, nrem_delta_waves);
%% Delta-spindle sequences
% Define spindle center time points
spindle_center    = round(sData.ephysdata2.NREMAbsSpindleIdx);

% Remove detected slow waves that are inside sleep spindles
slow_wave_start_stop = [slow_wave(:,1), slow_wave(:,4)];
log_idx = [];
for i = 1:length(spindle_center)
   log_idx = horzcat(log_idx, slow_wave_start_stop(:,1) > ECoG_spindle(i,1) & ...
       slow_wave_start_stop(:,2) < ECoG_spindle(i,2) );
end
slow_wave_removed = ~logical(sum(log_idx,2));

% Find delta troughs
[delta_troughs, sortId] = sort(slow_wave(slow_wave_removed,3), 'ascend');

% For each delta wave trough, define a time window from + 100ms to 1.3s in
% which to look for any spindle centers. 
delta_win_start   = srate*0.1;
delta_win_end     = srate*1.3;
delta_spindle_win = [delta_troughs+delta_win_start, delta_troughs+delta_win_end];

slow_wave = slow_wave(slow_wave_removed,: );
slow_wave = slow_wave(sortId,:); 
% delta_troughs = delta_troughs(slow_wave_removed(:,1));
% Preallocate
delta_spindle_idx_mat = zeros(  length(spindle_center), size(delta_spindle_win,1));

% Loop over time windows, define the specific time window for each delta
% wave trough, and check if any spindle center occurs within it. 
for i = 1:length(delta_spindle_win)
    time_win                   = delta_spindle_win(i,1):delta_spindle_win(i,2);
    delta_spindle_idx_mat(:,i) = ismember(spindle_center, time_win);
end

delta_spindle_spindle_idx  =  sum(delta_spindle_idx_mat,2);
enumerate_delta            = 1:length(slow_wave);
enumerate_spindles         = 1:length(ECoG_spindle);


%% Test
% delta_spindle_delta_idx   = logical( sum(delta_spindle_idx_mat,1));
% delta_trough_times = delta_troughs(delta_spindle_delta_idx);
% temp = delta_spindle_spindle_idx > 0;
% spindle_times = spindle_center(temp);
% close_delta_idx = zeros(1, length(delta_trough_times));
% for i = 1:length(close_delta_idx)
%     [~, close_delta_idx(i)] = min( abs(delta_trough_times(i)- spindle_times));
% end


% Because there can be multiple delta waves before a spindle, keep only
% first delta wave to avoid complicated indexing. 

% Find rows where there are multiple delta waves per spindle
row_idx  = delta_spindle_spindle_idx > 1;
row_nr   = enumerate_spindles(row_idx);
get_rows = delta_spindle_idx_mat(row_idx,:);

% Preallocate
delta_waves_to_keep = zeros( size(get_rows,1),1);

for i = 1:size(get_rows,1)
    delta_waves_to_keep(i)                                  = find(get_rows(i,:),1,'last');
    delta_spindle_idx_mat(row_nr(i),:)                      = 0;
    delta_spindle_idx_mat(row_nr(i),delta_waves_to_keep(i)) = 1;
end

% Find the indicies of delta-spindle delta waves and sleep spindles
delta_spindle_spindle_idx =  logical( sum(delta_spindle_idx_mat,2));
delta_spindle_delta_idx   = logical( sum(delta_spindle_idx_mat,1));

delta_spindle_seq_spindles = enumerate_spindles(delta_spindle_spindle_idx);
delta_spindle_seq_delta    = enumerate_delta(delta_spindle_delta_idx);
%% SWR-delta sequences
% Find SWR peaks
swr_idx = sData.ephysdata.absRipIdx';

% For each SWR peak, define a time window from + 50ms to 250ms in
% which to look for any delta troughs.
swr_win_start = srate*0.05;
swr_win_end   = srate*0.25;
swr_delta_win = [swr_idx+swr_win_start, swr_idx+swr_win_end];

% Preallocate
swr_delta_idx_mat = zeros( length(delta_troughs), size(swr_delta_win,1) );

% Loop over time windows and check if any delta wave falls within it. 
for i = 1:length(swr_delta_win)
    time_win               = swr_delta_win(i,1):swr_delta_win(i,2);
    swr_delta_idx_mat(:,i) = ismember(delta_troughs, time_win);
end

% Find the indicies of delta-spindle delta waves and sleep spindles
swr_delta_delta_idx = logical( sum(swr_delta_idx_mat,2));
swr_delta_swr_idx   = logical( sum(swr_delta_idx_mat,1));

enumerate_swr       = 1:length(swr_idx);

swr_delta_seq_delta = enumerate_delta(swr_delta_delta_idx);
swr_delta_seq_swr   = enumerate_swr(swr_delta_swr_idx);

%% Correct for SWR-delta mismatch
% Note that there can be a mismatch between nr of swrs and nr of delta
% waves in the SWR-delta sequences (Typically more SWRs than delta waves).
% Because only one SWR,delta, or spindle event is needed in the
% SWR-delta or SWR-delta-spindle sequence, if there are two SWRs linked to one delta
% wave, keep only the one that is closest to delta wave trough. This also
% makes indexing simpler. 

% First find SWR time points of those SWRs coupled to delta waves
swr_times = swr_idx(swr_delta_seq_swr);

% Then find delta wave time points of those waves coupled to SWRs. 
delta_times = delta_troughs(swr_delta_delta_idx);

% Loop over delta wave times and check which coupled SWR is closest. That
% SWR is defined as the delta wave associated SWR. 
close_swr_idx = zeros(1, length(delta_times));
for i = 1:length(delta_times)
    [~, close_swr_idx(i)] = min( abs(delta_times(i)- swr_times));
end

swr_delta_seq_swr = swr_delta_seq_swr(close_swr_idx);

%% SWR-delta-spindle sequences
% SWR-delta-spindle sequences are defined as the conjuction of
% delta-spindle and SWR-delta sequences.

% Find overlapping delta waves (i.e. waves that are both coupled to sleep 
% spindles AND to SWRs. 
delta_overlap_spindle_swr = ismember(delta_spindle_seq_delta, swr_delta_seq_delta);
delta_overlap_swr_idx     = ismember(swr_delta_seq_delta, delta_spindle_seq_delta );
swr_to_include            = swr_delta_seq_swr(delta_overlap_swr_idx);

%% Delta-SWR sequences
delta_swr_win_start = srate*0.05;
delta_swr_win_end   = srate*0.4;
delta_swr_win       = [delta_troughs+delta_swr_win_start, delta_troughs+delta_swr_win_end];
delta_swr_idx       = zeros(length(swr_idx), size(delta_swr_win,1));

for i = 1:length(delta_swr_win)
    time_win           = delta_swr_win(i,1):delta_swr_win(i,2);
    delta_swr_idx(:,i) = ismember(swr_idx, time_win);
end
delta_swr_idx = logical( sum(delta_swr_idx,2));

%% SWR-spindle sequences
swr_spindle_win_start = srate*0.05;
swr_spindle_win_end   = srate*1;
swr_spindle_win       = [swr_idx+swr_spindle_win_start, swr_idx+swr_spindle_win_end];
swr_spindle_idx_mat       = zeros(length(spindle_center), size(swr_spindle_win,1));

for i = 1:length(swr_spindle_win)
    time_win           = swr_spindle_win(i,1):swr_spindle_win(i,2);
    swr_spindle_idx_mat(:,i) = ismember(ECoG_spindle(:,1), time_win);
end
swr_spindle_spindle_idx = logical( sum(swr_spindle_idx_mat,2));
swr_spindle_swr_idx     = logical( sum(swr_spindle_idx_mat,1));
%% Output
% Delta-spindle events
delta_spindle_seq_spindles;
delta_spindle_seq_delta;

% SWR-delta events
swr_delta_seq_delta;
swr_delta_seq_swr;

% SWR-delta-spindle events
swr_delta_spindle_seq_swr     = swr_to_include;
swr_delta_spindle_seq_delta   = delta_spindle_seq_delta(delta_overlap_spindle_swr);
swr_delta_spindle_seq_spindle = delta_spindle_seq_spindles(delta_overlap_spindle_swr);

% Event rate 
n_delta_spindle_events     = length(delta_spindle_seq_spindles);
n_swr_delta_events         = length(swr_delta_seq_delta);
n_swr_delta_spindle_events = length(swr_to_include);
n_swr_spindle_events       = length(enumerate_spindles(swr_spindle_spindle_idx));

total_nrem  = sData.totalNREMmin;

nrem_n_delta_spindle     = n_delta_spindle_events/total_nrem;
nrem_n_swr_delta         = n_swr_delta_events/total_nrem;
nrem_n_swr_delta_spindle = n_swr_delta_spindle_events/total_nrem;
nrem_n_swr_spindle       = n_swr_spindle_events/total_nrem;

output = {nrem_n_delta_spindle, nrem_n_swr_delta, nrem_n_swr_delta_spindle,...
    nrem_n_swr_spindle};
%% Plot results

if length(varargin) > 1
    % Delta-spindle sequences    
%     if seq_to_plot == 1
        figure,
        plot(time,ecog_delta)
        hold on
        plot(time,ecog*.4-.2)
%         for i = enumerate_spindles(swr_spindle_spindle_idx)
%             z = [ECoG_spindle(i,1) ECoG_spindle(i,1) ECoG_spindle(i,2) ECoG_spindle(i,2) ]/2500;
%             v = [-.5 .5 .5 -.5];
%             patch(z, v, 'red', 'edgecolor', 'none', 'FaceAlpha', .3);
%         end
%         for i = enumerate_swr(swr_spindle_swr_idx)
%              z = [swr_start_stop(i,1) swr_start_stop(i,1) swr_start_stop(i,2) swr_start_stop(i,2) ]/2500;
%             v = [-.5 .5 .5 -.5];
%             patch(z, v, 'red', 'edgecolor', 'none', 'FaceAlpha', .3);
%         end
        for i = delta_spindle_seq_spindles
            z = [ECoG_spindle(i,1) ECoG_spindle(i,1) ECoG_spindle(i,2) ECoG_spindle(i,2) ]/2500;
            v = [-.5 .5 .5 -.5];
            patch(z, v, 'red', 'edgecolor', 'none', 'FaceAlpha', .3);
        end

        for i = delta_spindle_seq_delta
            x = [slow_wave(i,1) slow_wave(i,1) slow_wave(i,3) slow_wave(i,3)]/2500;
            y = [-.5 .5 .5 -.5];
            patch(x, y, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .3);
        end
        set(gca,'ylim',[-.5 .5])

    % SWR-delta sequences
%     elseif seq_to_plot == 2
%         figure,
%         plot(time,signal_filt)
%         hold on
%         plot(time,signal*.4-.1)
%         plot(time,hpc_lfp *.3 +.3)

        for i = swr_delta_seq_swr
            z = [swr_start_stop(i,1) swr_start_stop(i,1) swr_start_stop(i,2) swr_start_stop(i,2) ]/2500;
            v = [-.5 .5 .5 -.5];
            patch(z, v, 'black', 'edgecolor', 'none', 'FaceAlpha', .5);
        end

        for i = swr_delta_seq_delta
            x = [slow_wave(i,1) slow_wave(i,1) slow_wave(i,3) slow_wave(i,3)]/2500;
            y = [-.5 .5 .5 -.5];
            patch(x, y, 'cyan', 'edgecolor', 'none', 'FaceAlpha', .3);
        end
        set(gca,'ylim',[-.5 .5])

    % SWR-delta-spindle sequences
%     elseif seq_to_plot == 3
%         figure,
%         plot(time,signal_filt)
%         hold on
%         plot(time,signal*.4-.2)
        plot(time,hpc_lfp *.3 +.3)
        
        for i = swr_delta_spindle_seq_spindle
            z = [ECoG_spindle(i,1) ECoG_spindle(i,1) ECoG_spindle(i,2) ECoG_spindle(i,2) ]/2500;
            v = [-.5 .5 .5 -.5];
            patch(z, v, [0.6350 0.0780 0.1840], 'edgecolor', 'none', 'FaceAlpha', .3);
        end
        for i = swr_delta_spindle_seq_swr
            z = [swr_start_stop(i,1) swr_start_stop(i,1) swr_start_stop(i,2) swr_start_stop(i,2) ]/2500;
            v = [-.5 .5 .5 -.5];
            patch(z, v, [0.4660 0.6740 0.1880], 'edgecolor', 'none', 'FaceAlpha', .3);
        end
        for i = swr_delta_spindle_seq_delta
            x = [slow_wave(i,1) slow_wave(i,1) slow_wave(i,3) slow_wave(i,3)]/2500;
            y = [-.5 .5 .5 -.5];
            patch(x, y, [0 0.4470 0.7410], 'edgecolor', 'none', 'FaceAlpha', .3);
        end
        set(gca,'ylim',[-.5 .5])

%     end
end


