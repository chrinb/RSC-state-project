function plot_single_roi_spindle_activity(varargin)

%Written by Christoffer Berge | Vervaeke Lab

% Plot the activity of a single ROI across all spindles and its mean over all
% spindles. 
sData = varargin{1,1};
roinr = varargin{1,2};

nr_of_seconds   = 3;

% Add roisignals folder to path to get roi array
try
    folder_name = dir;
    addpath( [folder_name(1).folder '\roisignals'])
catch
end

% If a third input argument called 'axon' is given, merge highly correlated
% clusters
if nargin > 2
    % Find indicies of principal cells
    [pc_rois, in_rois] = remove_cells;
    if contains(varargin{1,3},'axon')
%             ROImtx                     = sData.imdata.roiSignals(2).newdff(pc_rois,:);
%             [dff,deconv, roiClustIDs, n_unchanged_rois] = hierClust_axons(ROImtx,1, sData, pc_rois);
        dff    = sData.imdata.roiSignals(2).mergedAxonsDff;
        deconv = sData.imdata.roiSignals(2).mergedAxonsDec; 
    elseif contains(varargin{1,3},'pc')
        dff    = sData.imdata.roiSignals(2).newdff(pc_rois,:);
        deconv = sData.imdata.roiSignals(2).ciaDeconvolved(pc_rois,:);
    elseif contains(varargin{1,3},'in')
        dff    = sData.imdata.roiSignals(2).newdff(in_rois,:);
        deconv = sData.imdata.roiSignals(2).ciaDeconvolved(in_rois,:);

    end
else
    dff   = sData.imdata.roiSignals(2).newdff;
    deconv = sData.imdata.roiSignals(2).ciaDeconvolved;
end


win_length      = (nr_of_seconds*31)*2+1;
time            = linspace(-nr_of_seconds,nr_of_seconds,win_length);


% Select the activity for a given ROI
% ROI_npil_signal   = npil(roinr, :);
roi_dff           = dff(roinr,:);
roi_dff_zscore    = zscore(roi_dff); 
roi_deconv        = deconv(roinr,:);
roi_deconv_zscore = zscore(roi_deconv);

% Choose spindle 
% prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
% spindle_band_select = input(prompt);
% 
% if spindle_band_select == 1
%     spindle_select = [];
% elseif spindle_band_select == 2
%     spindle_select = '1016';
% end
frames             = sData.daqdata.frame_onset_reference_frame;
spindle_select     = [];
spin_start_end_str = strcat('NREMspindleStartEnd', spindle_select);
spindleIdxStart    = frames(round(sData.ephysdata2.(spin_start_end_str)(:,1)));

% Preallocate
[dFF_spindle_activity,deconv_spindle_activity,dFF_spindle_activity_zscore,deconv_spindle_activity_zscore]...
    = deal(zeros(length(spindleIdxStart), win_length));

%% Specify whether to threshold deconvolved dF/F or not
   
% Loop over sleep spindles and extract 2P segments
for spindle_nr = 1:length(spindleIdxStart)
    spindle_window_start = spindleIdxStart(spindle_nr) - (nr_of_seconds*31); 
    spindle_window_end   = spindleIdxStart(spindle_nr) + (nr_of_seconds*31); 

    % skip spindles at the beginning or end with a window shorter than
    % length specified in seconds by user above
    if spindle_window_start > 1 && spindle_window_end < length(dff)
        dFF_spindle_activity(spindle_nr, :) = ...
            roi_dff(spindle_window_start:spindle_window_end); 

        dFF_spindle_activity_zscore(spindle_nr, :) = ...
            roi_dff_zscore(spindle_window_start:spindle_window_end); 

        deconv_spindle_activity(spindle_nr, :) = ...
            roi_deconv(spindle_window_start:spindle_window_end);

        deconv_spindle_activity_zscore(spindle_nr, :) = ...
            roi_deconv_zscore(spindle_window_start:spindle_window_end);
    end
end
% 
% if threshold == 1 
%     [threshold_percentile, ~] = threshold_deconvolved(deconv, roi_deconv, roinr);
%     ROI_idx = reshape(deconv_spindle_activity, 1,[]);
%     ROI_idx(ROI_idx < threshold_percentile(roinr)) = 0;
%     deconv_spindle_activity = reshape(ROI_idx, size(deconv_spindle_activity));
% end
deconv_spindle_activity        = smoothdata(deconv_spindle_activity,2, 'gaussian', 5);
deconv_spindle_activity_zscore = smoothdata(deconv_spindle_activity_zscore,2, 'gaussian', 5);
%% Plot results

c_lim_dff = [-3 3];
c_lim_dec = [0 2];
sessionID = sData.sessionInfo.sessionID;

x1 = [time(1) time(end)];
y1 = [1 length(spindleIdxStart) ];

figure,
sgtitle([ sessionID, ', ROI # ' num2str(roinr)], 'Interpreter', 'none')

subplot(241)
imagesc(x1,y1, dFF_spindle_activity)
xlabel('Time from spindle onset (sec)')
ylabel('Spindle #')
title('dF/F')

subplot(242)
imagesc(x1, y1, dFF_spindle_activity_zscore)
xlabel('Time from spindle onset (sec)')
ylabel('Spindle #')
title('z-score dF/F')

subplot(243)
imagesc(x1, y1, deconv_spindle_activity)
xlabel('Time from spindle onset (sec)')
ylabel('Spindle #')
title('deconvolved dF/F')


subplot(244)
imagesc(x1, y1, deconv_spindle_activity_zscore)
xlabel('Time from spindle onset (sec)')
ylabel('Spindle #')
title('z-score deconvolved dF/F')


subplot(245),
shadedErrorBar(time, nanmean(dFF_spindle_activity), ...
    std(dFF_spindle_activity)/sqrt(numel(dFF_spindle_activity(:, 1))) ,'lineprops', 'b');
xlabel('Time from spindle onset (sec)')
ylabel('Mean dF/F')
set(gca, 'xlim',[time(1) time(end)])

subplot(246),
shadedErrorBar(time, nanmean(dFF_spindle_activity_zscore), ...
    std(dFF_spindle_activity_zscore)/sqrt(numel(dFF_spindle_activity_zscore(:, 1))) ,'lineprops', 'b');
xlabel('Time from spindle onset (sec)')
ylabel('Mean z-score dF/F')
set(gca, 'xlim',[time(1) time(end)])

subplot(247),
shadedErrorBar(time, nanmean(deconv_spindle_activity), ...
    std(deconv_spindle_activity)/sqrt(numel(deconv_spindle_activity(:, 1))) ,'lineprops', 'b');
xlabel('Time from spindle onset (sec)')
ylabel('Mean deconvolved dF/F')
set(gca, 'xlim',[time(1) time(end)])

subplot(248),
shadedErrorBar(time, nanmean(deconv_spindle_activity_zscore), ...
    std(deconv_spindle_activity_zscore)/sqrt(numel(deconv_spindle_activity_zscore(:, 1))) ,'lineprops', 'b');
xlabel('Time from spindle onset (sec)')
ylabel('Mean z-score deconvolved dF/F')
set(gca, 'xlim',[time(1) time(end)])

x = roinr;
txt = ['ROI #', num2str(x)];
sgtitle(txt)