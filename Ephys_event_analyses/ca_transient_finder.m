function sData = ca_transient_finder(sData, params)

% Written by Christoffer Berge | Vervaeke lab

% Function that computes significant Ca2+ transients based on method
% described in Pettit et al. (2022).



%% Get exctitatory and inhibitory indices
[pc_rois, in_rois] = remove_cells_longitudinal(sData, params);

% Select data
switch params.cell_type
    case 'axon'
    dff    = sData.imdata.roiSignals(2).mergedAxonsDffFilt;
    txt    = 'mergedAxons';
    % case {'pc', 'in'}
    % dff = sData.imdata.roiSignals(2).newdff;
    % txt  = 'all';

    case 'in'
    dff = sData.imdata.roiSignals(2).newdff(in_rois,:);
    txt = 'in';
    case 'pc'
    dff = sData.imdata.roiSignals(2).newdff(pc_rois,:);
    txt = 'pc';
%     case 'all'
%     dff = sData.imdata.roiSignals(2).newdff;
end

% Standardize data (for each trace, subtract median and divide by standard
% deviation)

% For excitatory cells, the median provides a baseline DF/F at around zero
% and with little overestimation of significant transients. However, it
% underestimates a lot of transients for inhibitory cells (baseline is 
% typically above 0). 20th percentile might therefore be a better choice for INs

dff_stand = (dff - median(dff,2)) ./ std(dff,[],2); 
% dff_stand = (dff - prctile(dff,20,2)) ./ std(dff,[],2);

% Create threshold levels (in units of standard deviation)
std_thresholds = 1:.2:4;

%% Compute significant transients

% Preallocate
[sig_transients_logmat, sig_transients] = deal( zeros( size(dff)));

% Loop over ROIs
parfor roi_nr = 1:size(dff,1)
    tic;

    roi_signal      = dff_stand(roi_nr,:);
    roi_signal_orig = dff(roi_nr,:);

    % Loop over threshold levels
    threshold_vec = [];
    frame_vec     = [];
    for threshold_level = std_thresholds
        
        % Find positive samples exceeding current threshold (putative transients) and the nr of
        % frames of each (width)
        logical_pos_peaks = double(roi_signal > threshold_level);

        % Correct for ongoing transients occuring at imaging start
        if logical_pos_peaks(1) == true
            logical_pos_peaks(1) = false;
        end

        % Find the location and widths (nr of frames) of peaks exceeding current threshold
        [~, locs, width_pos, ~] = findpeaks(logical_pos_peaks, 'MinPeakHeight', 0, 'WidthReference', 'halfheight');
        
        % Find the number of negative-going transients (<-threshold_level)
        % at current threshold and widths (nr of frames) of each
        logical_neg_peaks    = double(roi_signal < -threshold_level);
        [~, ~, width_neg, ~] = findpeaks( logical_neg_peaks, 'MinPeakHeight', 0, 'WidthReference', 'halfheight'); 
        
        % Find unique nr of transient frames (n)
        unique_n = unique( width_pos);

        % Loop over putative transients of n frames in length
        new_vec   = [];
        for n_frames = 1:numel(unique( width_pos))
            
            % Select current transient n frame
            n_frame  = unique_n(n_frames);

            % Compute false positive rate
            neg            = sum(width_neg >= n_frame);
            pos            = sum(width_pos >= n_frame);
            false_pos_rate = neg/pos;
            
            % If false positive rate of transients of frame nr lenght = n
            % is less than 0.001 (0.1%), save the indicies of those
            % transients and logical vector
            if false_pos_rate < 0.001
                
                % Find indicies of significant transients and locations
                temp_idx  = width_pos == n_frame;
                ind       = locs(temp_idx);

                idx_vec      = zeros(1, length(roi_signal));
                idx_vec(ind) = 1;
                
                n_of_ones  = ones(1, n_frame-1);
                numel_skip = size(n_of_ones,2)+1;
                
                % Loop over significant transients time points (indicies) 
                % and insert logical ones for transient duration
                for transient_idx = 1:size(ind,2)
                    temp_vec = cat(2, idx_vec(1:ind(transient_idx)), n_of_ones,  idx_vec(ind(transient_idx)+numel_skip:end) );
                    new_vec  = [new_vec; temp_vec]; 
                end
               
                if ~( size(new_vec,1) == 1)
                    new_vec = sum(new_vec);
                end
                new_vec = logical(new_vec);
            end
            % Create matrix where each row is a logical vector containing
            % all transients of width n (n_frames) with a false positive rate > 0.001
            threshold_vec = [threshold_vec; new_vec];
        end

        % After computing all significant transients at current threshold
        % level, create new logical vector

        % If threshold_vec is empty, do nothing

        % If there is only one vector at current threshold level, convert
        % that to logical vector
        if size(threshold_vec,1) == 1
            threshold_vec = logical(threshold_vec);
        end
        % If more threshold_vec is matrix
        if size(threshold_vec,1) > 1
            threshold_vec = sum(threshold_vec);
        end
        
        % Store all logical vectors at current threshold level in new
        % matrix called frame_vec
        frame_vec = [frame_vec; threshold_vec];
    end

    frame_vec = sum(frame_vec);
    frame_vec = logical(frame_vec);
    
    [eventStartIdx, eventStopIdx ]  = findTransitions( frame_vec );
%     figure,
%     plot(frame_vec), hold on
%     set(gca, 'xlim',[1 200], 'ylim',[-1 2])

    % Merge transients separated by less than three frames 
    for n_transients = 1:size(eventStartIdx,2)-1
        if eventStartIdx(n_transients+1) - eventStopIdx(n_transients) <= 3
            frame_vec( eventStopIdx(n_transients):eventStartIdx(n_transients+1)) = 1;
        end
    end
    
    [eventStartIdx, eventStopIdx ]  = findTransitions( frame_vec );
%     figure, 
%     plot(frame_vec), hold on
%     set(gca, 'xlim',[1 200], 'ylim',[-1 2])
         
    % Remove transients less than two frames in duration
    for n_transients = 1:size(eventStartIdx,2)
        if length( eventStartIdx(n_transients):eventStopIdx(n_transients) ) == 1
            frame_vec( eventStartIdx(n_transients):eventStopIdx(n_transients)) = 0;
        end
%         plot(frame_vec, 'linew',1)
    end
    
    % Store logical indices of significant transients
    sig_transients_logmat(roi_nr,:) = frame_vec;
    
    % Set frames outside sig. transients to zero
    roi_sig_transients    = roi_signal_orig;
    roi_sig_transients(frame_vec == 0) = 0;
    
    % Store sig. transients in matrix
    sig_transients(roi_nr, :) = roi_sig_transients;

    % Set breakpoint here to plot the significant transients of current ROI

%     temp_plot_sig    = roi_signal;
%     temp_plot_nonsig = roi_signal;
%     temp_plot_sig(frame_vec == 0) = NaN;
% %     roi_signal_sig_transients(frame_vec == 0) = NaN;
%     figure(1), 
%     hold on
%     grayColor = [.6 .6 .6];
%     plot(temp_plot_nonsig, 'Color', grayColor)
%     plot(temp_plot_sig, 'r')
%     clf
    t = toc;
    fprintf('\n Finding significant transients for ROI # %d in %.1f seconds', roi_nr,t)

end

%% Store output
sData.analysis.transients.([params.cell_type, '_sig_transients_logmat']) = sig_transients_logmat;
sData.imdata.roiSignals(2).([txt, '_sig_transients'])                     = sig_transients;

%% Optional: plot some example significant traces

% n_rois_to_plot = 10;
% time_imaging   = linspace(1, size(dff_stand,2), size(dff_stand,2) )/31;
% 
% switch varargin{1,2}
%     case 'pc'
%     rois_to_plot = randsample(pc_rois, n_rois_to_plot);
%     case 'in'
%     rois_to_plot = randsample(in_rois, n_rois_to_plot);
% %     case 'axon'
% %     rois_to_plot = 
% end
% 
% figure,
% 
% 
% max_lim   = max( max(dff(rois_to_plot,:)));
% min_lim   = min( min(dff(rois_to_plot,:)));
% 
% grayColor = [.6 .6 .6];
% % timescale = dsearchn(time_imaging', [0 30]');
% for i = 1:n_rois_to_plot
% 
%     roi = rois_to_plot(i);
%     % Plot scale bars
%     if i == n_rois_to_plot
%         h(i) = subplot( n_rois_to_plot ,1,i);
%         plot(time_imaging, ones(1, size(dff,2)), 'Color', [1 1 1] ), hold on
%         plot([0 0], [1 0], 'LineWidth',3, 'Color',[0 0 0]), hold on
%         plot([0 30], [0 0], 'LineWidth', 3, 'Color',[0 0 0])
%         set(gca,'visible','off', 'ylim', [min_lim max_lim])
%     else
%     h(i) = subplot( n_rois_to_plot ,1,i);
%     plot(time_imaging, dff(roi,:), 'Color', grayColor)
%     hold on
%     plot(time_imaging, sData.analysis.transients.roi_signal_sig_transients(roi,:),'r')
% 
%     set(gca,'visible','off', 'xlim', [time_imaging(1) time_imaging(end)],...
%         'ylim', [min_lim max_lim])
%     end
% 
% end
% linkaxes(h, 'x')
% 
% 
