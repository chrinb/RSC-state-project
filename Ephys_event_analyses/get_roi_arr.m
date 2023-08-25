function sData  = get_roi_arr(sessionObject)

% Written by Christoffer Berge | Vervaeke lab

% Function that adds roi array to sData struct

sessionObject;
% Add roisignals folder to path to get roi array
try
    folder_name = dir;
    addpath( [folder_name(1).folder '\roisignals'])
catch
end

folder_name = dir;
temp_var = folder_name(1).folder;
roiFolder = [temp_var '\roisignals'];
filesInFolder = dir(roiFolder);

% Loop over files in folder
for x = 1:length(filesInFolder) 
    if ~isempty(strfind(filesInFolder(x).name,'rois')) == 1
        try 
            load(filesInFolder(x).name);
        catch
        end
    end
end

sData.imdata.roiArr = roi_arr;