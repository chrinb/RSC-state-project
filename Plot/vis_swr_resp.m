function interesting_swr_idx = vis_swr_resp(sData)

% Function that allows user to visualize the activity of ROI's across
% individual SWRs

nr_of_seconds = 3;
win_length    = (nr_of_seconds*31)*2+1;
time          = linspace(-nr_of_seconds,nr_of_seconds,win_length);
frames        = sData.daqdata.frame_onset_reference_frame;
%% Select 2P signal
signal_dff         = sData.imdata.roiSignals(2).newdff;
signal_deconvolved = sData.imdata.roiSignals(2).ciaDeconvolved;

%% Select SWRs for analysis
prompt = sprintf('Select SWRs to analyze: 1 = All | 2 = Awake | 3 Spindle-Uncoupled | 4 Spindle-Coupled ');
swrs_to_analyze = input(prompt);

if swrs_to_analyze == 1 
    swr_idx = sData.ephysdata.absRipIdx;
elseif swrs_to_analyze == 2
    swr_idx = sData.ephysdata.awake_swr;
elseif swrs_to_analyze == 3 
    swr_idx = sData.ephysdata.NREM_spindle_uncoupled_swr;
elseif swrs_to_analyze == 4 
    swr_idx = sData.ephysdata.spindle_coupled_swr;
end

prompt = sprintf('Remove locomotion SWR? (1 = yes) ');
rem_run_swr = input(prompt);
    
prompt = sprintf('Remove temporally close SWR? (1 = yes) ');
rem_cluster_swr = input(prompt);
    
if ~isempty(rem_run_swr) && isempty(rem_cluster_swr)
    [swr_idx,~] = riprun2(sData, swr_idx); 
elseif isempty(rem_run_swr) && ~isempty(rem_cluster_swr)
    swr_idx = RemoveRip(swr_idx);
elseif ~isempty(rem_run_swr) && ~isempty(rem_cluster_swr)
   [swr_idx,~] = riprun2(sData, swr_idx); 
   swr_idx = RemoveRip(swr_idx);
end

% time_vec = (-3:/
% Convert e-phys SWR times to 2P frame rate times
swr_idx = frames(round(swr_idx));
t = 1;
x1 = [time(1) time(end)];
y1 = [1 size(signal_dff,1) ];
swr_nr = 1;
while swr_nr <= length(swr_idx)
    swr_window_start = swr_idx(swr_nr) - (nr_of_seconds*31); 
    swr_window_end   = swr_idx(swr_nr) + (nr_of_seconds*31); 
    if swr_window_start > 1 && swr_window_end < length(signal_dff)
        figure,
        subplot(211)
        imagesc(x1, y1, signal_dff(:,swr_window_start:swr_window_end));
        colormap jet
        colorbar
        caxis([0 1])

        subplot(212)
        imagesc(x1, y1, signal_deconvolved(:,swr_window_start:swr_window_end));
        colorbar
        
        txt = ['SWR #', num2str(swr_nr)];
        sgtitle(txt)
        
        prompt = sprintf('Store SWR index?');
        x = input(prompt,'s');
    
        if strcmp(x,'y')
            interesting_swr_idx(t) = swr_nr;
            t = t+1;
            close
            swr_nr = swr_nr + 1;
        elseif strcmp(x,'b')
            close
            swr_nr = swr_nr - 1;
        else
            close
            swr_nr = swr_nr + 1;
        end
    else
        swr_nr = swr_nr + 1;
    end
end

interesting_swr_idx;
    