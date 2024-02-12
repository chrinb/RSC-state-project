
srate = find_imaging_framerate(sData);

% Get event times for 1 neuron
dec = sData.imdata.roiSignals(2).ciaDeconvolved;
dff = sData.imdata.roiSignals(2).newdff;
% Get SWR events and convert to seconds
swr = sData.ephysdata.absRipIdx;
swr_frame = sData.ephysdata.frameRipIdx;
swrS = swr.*(1/2500);
swrF = swr_frame.*(1/srate);

temp1 = zeros(1, size(dec,2));

temp1(swr_frame ) = 1;

time_vec = (0:size(dec,2))/srate;

sData.daqdata.frame_onset_reference_frame(end)
signal_length = 1:size(dec,2);
timeS = signal_length/srate;

%% Run Zeta test over all cells

% Try moving event with 0.5 s
dblBaselineDuration = 3;
matEventTimesWithPrecedingBaseline = swrS' - dblBaselineDuration;
rng(1,'twister'); %to make the output deterministic
    
% neuron_nr = 1:size(dec,1)
% for roi_nr = 1:size(dec,1)
    roi_nr = 1;
    neuron1 = dec(roi_nr,:);
    
    dec1s_idx = neuron1 > 0;
    vecSpikeTimes = timeS(dec1s_idx)';

    intPlot = 0;
    dblUseMaxDur = 3; %minimum of trial-to-trial durations

    % hTic=tic;
    try
        % [dblZetaP_pb, sZETA_pb, sRate_pb, sLatencies_pb]  = zetatest(vecSpikeTimes,matEventTimesWithPrecedingBaseline,[],[],3,[],[],[],[]);
       
        dblZetaP_pb = zetatest(vecSpikeTimes,swrS,dblUseMaxDur,[],3);
        
        dblBaselineDurationMs = dblBaselineDuration*1000;
    
        drawnow;hFig = gcf;
        for intPlotNr=1:numel(hFig.Children)
	        %adjust x-ticks
	        if contains(hFig.Children(intPlotNr).XLabel.String,'Time ')
		        set(hFig.Children(intPlotNr),'xticklabel',cellfun(@(x) num2str(str2double(x)-dblBaselineDuration),get(hFig.Children(intPlotNr),'xticklabel'),'UniformOutput',false));
	        end
            %adjust timings in title
            strTitle = hFig.Children(intPlotNr).Title.String;
            [vecStart,vecStop]=regexp(strTitle,'[=].*?[m][s]');
            for intEntry=1:numel(vecStart)
                strOldNumber=hFig.Children(intPlotNr).Title.String((vecStart(intEntry)+1):(vecStop(intEntry)-2));
                strTitle = strrep(strTitle,strcat('=',strOldNumber,'ms'),strcat('=',num2str(str2double(strOldNumber)-dblBaselineDurationMs),'ms'));
            end
            hFig.Children(intPlotNr).Title.String = strTitle;
        end
        drawnow;
        
        %here we adjust the times in the variables that getZeta returns
        sLatencies_pb.Onset = sLatencies_pb.Onset - dblBaselineDuration;
        sLatencies_pb.Peak = sLatencies_pb.Peak - dblBaselineDuration;
        sLatencies_pb.ZETA = sLatencies_pb.ZETA - dblBaselineDuration;
        sLatencies_pb.ZETA_InvSign = sLatencies_pb.ZETA_InvSign - dblBaselineDuration;
        sZETA_pb.vecSpikeT = sZETA_pb.vecSpikeT - dblBaselineDuration;
        sZETA_pb.vecLatencies = sZETA_pb.vecLatencies - dblBaselineDuration;
        sZETA_pb.dblZetaT = sZETA_pb.dblZetaT - dblBaselineDuration;
        sZETA_pb.dblZetaT_InvSign = sZETA_pb.dblZetaT_InvSign - dblBaselineDuration;
        sRate_pb.vecT = sRate_pb.vecT - dblBaselineDuration;
        sRate_pb.dblPeakTime = sRate_pb.dblPeakTime - dblBaselineDuration;
        sRate_pb.dblOnset = sRate_pb.dblOnset - dblBaselineDuration;
    catch
        dblZetaP_pb = NaN;
    end
    zeta_output(roi_nr)= dblZetaP_pb;

   
    clear dblZetaP_pb sZETA_pb sRate_pb sLatencies_pb
    % dblElapsedTime=toc(hTic);
    % fprintf('\nDefault parameters (elapsed time: %.2f s):\nzeta-test p-value: %f\n',dblElapsedTime,dblZetaP_pb)
% end

[zsort, zsort_idx] = sortrows(zeta_output');