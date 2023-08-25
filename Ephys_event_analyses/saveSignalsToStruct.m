function [sData] = saveSignalsToStruct(sData)

% Written by Christoffer Berge | Vervaeke lab

% Function that saves extracted imaging data to sData struct
% after using "save signals" function in roimanager

% Add roisignals folder to path to get roi array
try
    folder_name = dir;
    addpath( [folder_name(1).folder '\roisignals'])
catch
end

% fileFolder = uigetdir('D:\Data');
folder_name = dir;
temp_var = folder_name(1).folder;
roiFolder = [temp_var '\roisignals'];
filesInFolder = dir(roiFolder);

for x = 1:length(filesInFolder) 
    if ~isempty(strfind(filesInFolder(x).name,'rois')) == 1
        try 
            load(filesInFolder(x).name);
        catch
        end
    end

    if ~isempty(strfind(filesInFolder(x).name,'signals')) == 1 % if the file is a signals.mat file
        try 
            load(filesInFolder(x).name)
        catch
            warning('Add "roisignals" to path')
        end          
    end    
end


% assign the variables into their appropriate slot in the sData struct
try sData.imdata.roiSignals(2).ciaDeconvolved = ciaDeconv; catch;end
try sData.imdata.roiSignals(2).ciaDenoised = ciaDenois; catch; end
sData.imdata.roiSignals(2).newdff = dff;
try sData.imdata.roiSignals(2).npilMediF = npilMediF; catch; end
try sData.imdata.roiSignals(2).roisMeanF = roisMeanF; catch; end
sData.imdata.roiSignals(2).roisMeanFRaw = roisMeanFRaw;

try sData.imdata.roi_arr = roi_arr; catch; end

% folder_name = dir(fileFolder);
folder_cell = strfind({folder_name.name}, '.mat');
for i = 1:length(folder_cell)
    if folder_cell{1,i} > 0
    idx = i;
    end
end

if length(idx) > 1
    error('Idx is > 1, could not locate sData file name')
end

fileName = folder_name(idx).name;
clearvars -except sData fileName
save(fileName)
