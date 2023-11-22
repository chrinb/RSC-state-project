function plot_meanF_DFF_comparison(sData, sessionObjects, roi_nr, denoise)


% Load mean fluorescence, DF/F
roi_meanF_table = sessionObjects.loadData('RoiSignals_MeanF');
roi_meanF       = table2array(roi_meanF_table); 
roi_meanF       = roi_meanF';

run_dff   = sData.imdata.roiSignals(2).newdff;
dec       = sData.imdata.roiSignals(2).ciaDeconvolved;
transient = sData.imdata.roiSignals(2).all_sig_transients;


grid_classification = sData.imdata.ch2_grid_classficiation;

imaging_sampling_rate = find_imaging_framerate(sData);

time_vec = (0:size(run_dff,2)-1)./imaging_sampling_rate;


%% Plot
if grid_classification(roi_nr) == 0
    msgbox('This ROI is excluded via grid analysis')
else

    if denoise == 1
        meanF   = okada(roi_meanF(roi_nr,:),2);
        meanDFF = okada( run_dff(roi_nr,:),2);
    else
        meanF   = roi_meanF(roi_nr,:);
        meanDFF = run_dff(roi_nr,:);
    end

    figure, 
    sgtitle(['ROI # ', num2str(roi_nr)])
    h(1) = subplot(211);
    plot(time_vec, meanF, 'k');
    ax = gca;
    ax.XAxis.Visible = 'off';
    ylabel('Mean fluorescence', FontSize=16)


    h(2) = subplot(212); hold on
    plot(time_vec, meanDFF, 'k');
    plot(time_vec, dec(roi_nr,:)-.3, 'r');
    xlabel('Time (s)', FontSize=16)
    ylabel('DF/F', FontSize=16)
    
    % dff = run_dff(roi_nr,:)
    % dff( sData.analysis.transients.pc_sig_transients_logmat(1,:) == 0)
    % transient_log = sData.analysis.transients.pc_sig_transients_logmat;
    % 
    % transient_signal = transient(roi_nr,:);
    % transient_signal(transient_signal == 0) = NaN;
    % 
    % h(3) = subplot(313);
    % plot(time_vec, meanDFF, 'k'); hold on
    % plot(time_vec, transient_signal, 'r')

end

linkaxes(h, 'x');
set(gca, 'xlim', [time_vec(1) time_vec(end)])