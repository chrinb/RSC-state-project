function plot_place_cells(sData, params, sessionObject)

%{
Plot colormap of place cells (cells x time) and wheel position below
%}

subjectID = sessionObject.subjectID;

if strcmp(params.signal_type, 'dff')
    dff         = sData.imdata.roiSignals(2).newdff;
    place_cells = dff(sData.imdata.MaoPC_dff.PlaceCellActSortingOrder,:);

    if strcmp(subjectID, 'm6157') || strcmp(subjectID, 'm6159')
        place_cells = okada(place_cells, 2);
    end
elseif strcmp(params.signal_type, 'dec')
    dec         = sData.imdata.roiSignals(2).ciaDeconvolved;
    place_cells = dec(sData.imdata.MaoPC_dff.PlaceCellActSortingOrder,:);
    place_cells(place_cells>0) = 1;
end


srate = find_imaging_framerate(sData);
time_imaging  = linspace(1, size(place_cells,2), size(place_cells,2) )/srate;
y1 = [1 size(place_cells,1)];

%% Plot 
figure(3), clf

h(1) = subplot(3,2, 1:4);
imagesc(time_imaging, y1, place_cells)
h(1).XTickLabel = [];
% caxis([0 0.4])
if strcmp(params.signal_type, 'dec')
    colormap(flipud(gray));
end
ylabel('Place cell #')
set(gca, "FontSize", 12)
title(sData.sessionInfo.sessionID);

h(2) = subplot(3,2,5:6);
plot(sData.behavior.timeInSec  , sData.behavior.wheelPosDs(1:end-1), 'k','LineWidth',1)
ylabel('Postion (cm)')
xlabel('Time (s)')
set(gca, 'xlim',[time_imaging(1) time_imaging(end)], 'FontSize', 12)
linkaxes(h, 'x')