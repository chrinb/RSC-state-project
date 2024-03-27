%%%%% NORA'S ANALYSIS PIPLINE FOR ACTIVITY OF CELLS (IMAGED CA-TRANSIENTS AND BEHAVIOR OF THE ANIMAL IN A SPATIAL NAVIAGTION TASK
                            
%%%%% MOTION CORRECTION OF CA-SESSIONS AND DATA EXTRACTION

%%%%% CONVERSION AND ANALYSIS OF LABVIEW DATA FOR BEHAVIOR

%%% Create mouseInfo struct into C\MATLAB\MOUSEINFO folder (once for each mouse)
% Use: CreateMouseSheet_empty.mat (not a function), save the file as such:
% mouseInfo-XXXX.mat, where XXXX is the mouse ID (e.g. 8001)

%%% SET:
% mouseFolder = 'C:\Users\chrinb\Desktop\PhD\Code\MOUSEINFO'; % save here all mouseInfo files
current_session_dir    = dir;
% mouseFolder = current_session_dir(1).folder(1:18); % Previous code
mouseFolder = current_session_dir(1).folder;

% Load sData

% savePath    = 'C:\Users\chrinb\Desktop\PhD\Code\SAVE'; % a new folder will be generated with the session name in this folder
savePath           = current_session_dir(1).folder;
mouseWeightPercent = 80; % set mouse weight in percentage compared to original weight before water restriction
nBins              = 80; % number of bins you want to divide your wheel circumference.

%% ANALYSING BEHAVIOR DATA:
% Create sessionInfo and mouseInfo in sData:

[fileName,filePath,~] = uigetfile('*.tdms','',[savePath, '\labview_behavior'] );
% load([fileName(1:17), '.mat']);
sData.sessionInfo.mouseWeightPercent = mouseWeightPercent;
sData.sessionInfo.savePath           = [mouseFolder, '\Behavior'];

ID_pt1                      = savePath(27:40);
ID_pt2                      = savePath(51:54);
sData.sessionInfo.sessionID = [ID_pt1, ID_pt2];
% [sData,filePath,fileName] = SA_CreateSessionInfo(mouseFolder,savePath,mouseWeightPercent); % mouseInfo is needed

%%% Convert raw labview behavioral tdms data into sData.daqdata struct (lick, position, reward, frame signal)
sData = SA_loadTDMSdata(sData,filePath,fileName);

%%% Calculate behavior data (binning, plots)
sData = SA_CalcBehav(sData,nBins); 

%%% Comparing behavior of the first - second half of the session
LapsTested = 20; % how many laps take to analyse in the beginning and end
sData      = SA_BehavFirstSecondHalf(sData,LapsTested); 

% save(fullfile(sData.sessionInfo.savePath(1:54),...
%     strcat(sData.sessionInfo.sessionID(1:14), sData.sessionInfo.sessionID(25:28),'.mat')),'sData');
% save(fullfile(sData.sessionInfo.savePath(1:54),...
%     strcat(sData.sessionInfo.sessionID,'.mat')),'sData');
%close all

%% ANALYSING CALCIUM IMAGING DATA:
% you have to draw ROIs and extract the data in ROI manager:
% SET:
VelMin       = 0.1; % speed limit in cm/s. Below this speed, discard calcium data,set 0.1 for pyramidal cells
pmtGainGreen = 25;
waveLength   = 930; % nm
laserPower   = 100; % mW
fovCoord     = [0 0]; % AP, ML coordinated from center of window 
IsDeconv     = 0; % 0 : do not use deconvolved data, use only dFF for analysis. 1: use deconvolved data (already done and extracted by ROI manager), 2: do the deconvolution using the script.  
gol          = 5; % task name, lately I used gol=5 task. Landmark positions are written in the code 

% if ~isfield(sData.imdata, 'binned')
[sData, dff] = SA_calcCaData(sData, VelMin, IsDeconv, gol); %calcCaDataNori %calcCaDataNoriShortVIP
% end
clearvars nBins VelMin  laserPower waveLength fovCoord
close all

%% Place cell analysis based on Mao 2017, he used deconvolved data, and I adapted the code for dFF signal 

% Final params for dFF data (datatype=0): activityTreshold = 0.5; MinPlaceFieldSize = 2; MaxPlaceFieldSize = 100; InOutRatio = 2; ReliabiliyIndex = 0.3;
% Parameters originally used by Mao for deconvolved data (datatype=1):
% activityTreshold = 0.3; MinPlaceFieldSize = 15; MaxPlaceFieldSize = 120;
% InOutRatio = 3; ReliabiliyIndex = 0.34; binSize = 1.5 cm (I use 2 cm bins)

% dataype = 0 use dff data , datatype = 1 use deconvolved data

dff = sData.imdata.roiSignals(2).newdff(sData.imdata.roi_classification == 1, :);

datatype = 0; 
sData    = SA_placeCellMaosData(sData,datatype, dff); 

% detects cells which respond to both cues (e.g. both hot glue spikes or both velcro), 
% detects place cells with only one place field as well
sData = SA_LandmarkCellDetection(sData); 
 
% Save file to same path where other files can be found 
% save(fullfile(sData.sessionInfo.savePath(1:54),...
%     strcat(sData.sessionInfo.sessionID(1:14), sData.sessionInfo.sessionID(25:28),'.mat')),'sData');

save(fullfile(sData.sessionInfo.savePath(1:54),...
    strcat(sData.sessionInfo.sessionID,'.mat')),'sData');

close all;

