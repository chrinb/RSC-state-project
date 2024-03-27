rmaps = struct();

for k =1:length(data)
    
    rmaps(k,j).dff  = [];
    rmaps(k,j).deconv  = [];

    for j = 1:length(data(k).sessionIDs)
        sData  = mData(k,j).sData;        
        signal_dff = sData.imdata.roiSignals(2).dff;
        signal_deconv = sData.imdata.roiSignals(2).deconv;

        for f = 1:size(signal_dff,1)                
                signalMatr                  = helper.splitInTrialsPos(signal_dff, sData, f);                
                rmaps(k,j).dff(f,:)        =  nanmean(signalMatr,1);

                signalMatr             = helper.splitInTrialsPos(signal_deconv, sData, f);
                [signalMatr, ~]        = smoothdata(signalMatr, 2,'gaussian',8);
                
                rmaps(k,j).deconv(f,:) = nanmean(signalMatr,1);

        end
        
    end
end