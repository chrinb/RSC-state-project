function [sData,filePath,fileName] = SA_CreateSessionInfo(mouseFolder,savePath,mouseWeightPercent)

% msgbox('Choose TDMSdata file for behavior');

% [fileName,filePath,~] = uigetfile('*.tdms','','C:\' ); 
[fileName,filePath,~] = uigetfile('*.tdms','',[savePath, '\labview_behavior'] );


% parts = strsplit(fileName,'.'); % Clear".tdms" from Filename.
% sessionID = parts{1};

%%% SESSIONINFO

% Load sData
load(['m',fileName(1:17), '.mat']);
% These fields are already given by the "pipe_loadSessionDataS" function

% sData.sessionInfo.sessionID = sessionID;
% sData.sessionInfo.sessionNumber = str2double(sessionID(16:18));
% sData.sessionInfo.recordedData = {'2P'};
sData.sessionInfo.mouseWeightPercent = mouseWeightPercent;

%%% MOUSEINFO
% mouseFolderPath = [mouseFolder '\mouseInfo-' sessionID(2:5) '.mat'];

% These fields are already given by the "pipe_loadSessionDataS" function
% mouseFolderPath = [mouseFolder 'mouseInfo.mat'];
% 
% if isfile(mouseFolderPath)
% %     load([mouseFolder '\mouseInfo-' sessionID(2:5)]); % fileName might change
%     load(mouseFolderPath);
%     sData.mouseInfo = mouseInfo;
%     clear('mouseInfo');
% else
%     msgbox(['1) Make sure if "mouseFolder" variable is defined in the script and refers to the path where the mouseinfo files are stored.',char(10),char(10),... 
%         '2) Make sure if the mouse info file for this mouse is filled in manually and saved according to the mouse naming standards.'],'Mouseinfo data is not found.');
% end
 
%%% SAVING
% parts = strsplit(sessionID,'-'); % short version of sessionID
% fileID = strcat(parts{1},'-',parts{2},'-',parts{3});
% sData.sessionInfo.fileID = fileID;
% mkdir(savePath,fileID);
% sData.sessionInfo.savePath = strcat(savePath,'\',fileID);
% save(fullfile(sData.sessionInfo.savePath,strcat(fileID,'_sData.mat')),'sData');

end