function [BinLickMatrixExtended,MeanLickBinExtended,LickCmMatrixExtended,MeanLickCmExtended] =  SA_plotHeatBinLicks(sData) 

%PARAMETERS TO BE SET:
trial_nr         = sData.behavior.wheelLap;
bin_nr           = sData.behavior.meta.nBins;
EnterIntoBin     = sData.behavior.binning.enterIntoBinIndexExtended;
SampleSpentInBin = sData.behavior.binning.SampleSpentInBinExtendedBins;
bin_nr_plus      = size(EnterIntoBin,2);

% READ DATA FROM FILE
Temp    = 0; % insert an extra zero to the beginnning
Licking = vertcat(Temp,sData.behavior.lickDs); % check if there is lick within at least half the duration if the frame scan, if yes, consider lick (=1)

% set trial number for plotting, search last full data trial
for i = 1:trial_nr
    if any(isnan(EnterIntoBin(i,:))) || any(EnterIntoBin(i,:)==0)
        trial_nr_plot = i-1; 
        break
    else
        trial_nr_plot = i;
    end
end


BinLickMatrixExtended = NaN(trial_nr_plot, bin_nr_plus);
% calculate the number of licks during each bin (lick/cm)
for i = 1:trial_nr_plot  % rows are trials
    for j = 1:1:bin_nr_plus  % columns (distance bins)  
        LR = 0; % set lickrate to zero before counting in each bin
        for m = 1:1:SampleSpentInBin(i,j) % Note: if the recorging starts with lick, it will not detect that first lick
            if Licking(EnterIntoBin(i,j) + m - 2) == 0 && Licking(EnterIntoBin(i,j) + m-1) == 1   % search the first timepoint (sample) when licking starts, zero changes to one, this will be a lick event
                LR = LR + 1;  % if there is a lick event, increase lick rate with one
            end
        end
        BinLickMatrixExtended(i,j) = LR; % put the cumulative lick number in a bin into BinLick matrix
    end
end
MeanLickBinExtended = mean(BinLickMatrixExtended,1, 'omitnan');

LickCmMatrixExtended = BinLickMatrixExtended./sData.behavior.meta.binSize;
MeanLickCmExtended   = mean(LickCmMatrixExtended,1, 'omitnan');
%Cmax = ceil(max(MeanLickCm)*2);

%PLOT FIGURE
figure('Color','white'); 

x_axis = 1:sData.behavior.meta.binSize:(bin_nr*sData.behavior.meta.binSize);
y_axis = 1:trial_nr_plot;

imagesc(x_axis, y_axis, (LickCmMatrixExtended(1:trial_nr_plot,1:bin_nr))) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))

j = colorbar;
colormap(jet);
j.Label.String   = 'Licking (licks/cm)';
j.Label.FontSize = 11;
j.TickDirection  = 'out'; 
caxis([0 2]); 
hold on;
line([157 157],[0 trial_nr_plot],'Color','white');
hold on
line([163 163],[0 trial_nr_plot],'Color','white');

xlabel('Position on wheel (cm)');
ax = gca;
ax.TickDir = 'out';
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

ylabel('Trials');
if trial_nr_plot < 50
    yticklabels = 0:5:trial_nr_plot;
elseif trial_nr_plot >= 50 && trial_nr_plot < 200
    yticklabels = 0:10:trial_nr_plot;
else
    yticklabels = 0:20:trial_nr_plot;
end

% These lines sometimes result in axes labels that don't match the actual
% trial nrs...
% yticks = linspace(1, trial_nr_plot, numel(yticklabels));
% set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(sData.sessionInfo.sessionID);

end