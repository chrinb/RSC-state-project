function plotMultipleROIdFFNormAbsDist(plotdata,FrameRate,nROIOnFig,AbsDist,Lick,savePath,titl) % ROIsignals_dFF_Dombeck   %LightOn,
% PLOT multiple ROIs, max = 50

%plotdata = ROIsignals_dFF_Dombeck; % ROIsignals_raw % ROIsignals_raw_slow_removed %ROIsignals_dFF_Dombeck %ROIsignals_dFF_lowpass
% ROINuOnFig = 49; % how many ROIs to show on a fig


[ROINu,SampleNu] = size(plotdata);
startSample = 1;
stopSample = SampleNu-1;
FigNu = ceil(ROINu / nROIOnFig);
TimeMin = (1:SampleNu)/FrameRate/60; 

for m = 1:1:FigNu
    startROI = m*nROIOnFig - nROIOnFig + 1;
    stopROI = m*nROIOnFig;
    if stopROI > ROINu
       stopROI = ROINu; 
    end
    f1 = figure(2); 
    clf(f1,'reset')
    f1.Position = [100 100 1000 1000];
    f1.Color = [1 1 1];
    f1.Visible = 'on';
    % figure,
    plot(TimeMin(1,startSample:stopSample),Lick(startSample:stopSample)/1.2-1,'LineWidth', 1.5); hold on; 
    plot(TimeMin(1,startSample:stopSample),AbsDist(startSample:stopSample)/100,'LineWidth', 1.5); hold on; 
    %plot(TimeMin(1,startSample:stopSample),LightOn(startSample:stopSample)/1.2-2,'LineWidth', 1.5); hold on;
    for i = startROI:1:stopROI 
        ROIMax = max(plotdata(i,:)); % max Y value in ROI , data has to be scale to it
        plotROI = (plotdata(i,:))/ROIMax + i - startROI + 2 ;
        plot(TimeMin(1,startSample:stopSample),plotROI(1,startSample:stopSample),'LineWidth', 1.5); hold on; 
    end
    xlabel('Time (min)');
    ax = gca;
    ax.TickDir = 'out';
    ylabel(sprintf('ROI %d - %d (normalized dF/F)',startROI,stopROI));
    set(gca,'YTickLabel',[]);
    title(titl);
    fname = strcat(titl,sprintf('-ROI-%d-%d',startROI,stopROI));
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    close(2)
end

end