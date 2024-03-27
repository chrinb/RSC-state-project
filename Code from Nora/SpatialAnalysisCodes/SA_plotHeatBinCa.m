function SA_plotHeatBinCa(data,fileID,roi,ylab,BinSize,BinNu,TRNu,visibility) 

% parameters: x=160; y=CADATA.TRNu; data=CADATA.VeloLimInBin;
% data: matrix in which rows are trials, columns are bins, binned velocity values in cells
% ylabel : if Velo set it to 'Speed (cm/s)' , if dff set it to 'dF/F'
% Xaxis = 1:2:155,Yaxis = 1:TRNu

if roi == 0
    figtitle = fileID;
else
    figtitle = strcat(fileID,sprintf(' ROI #%d',roi));
end


Xaxis = 1:BinSize:BinSize*BinNu;

%% Plot
f1 = figure(3); 
f1.Color = 'white';
f1.Visible = visibility;

% subplot(3,2, 1:4)

imagesc(Xaxis, 1:TRNu, data);
c = colorbar;
colormap(viridis);
c.Label.String = ylab;
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
caxis([0 inf]);

axis([0 152 1 TRNu]);
xlabel('Position on Wheel (cm)');
% set(gca,'xticklabel',[])    
xticks([0,25,50,75,100,125,150]);
ylabel('Trials');
title(figtitle);

axis square