function metadata = loadImagingMetadataNoriMCFolder(List)
%loadImagingMetadata Load imaging metadata for specified session.
%   METADATA = loadImagingMetadata(sessionID) returns a struct with
%   imaging metadata from specified session.
%
%   METADATA contains following fields:
%       - dt            :   interframe interval
%       - xpixels       :   width of image in pixels
%       - ypixels       :   height of image in pixels
%       - objective     :   not necessary for now (not added)
%       - zoomFactor    :   zoomfactor during the recording
%       - pmtCh2        :   missing - can we implement??
%       - pmtCh3        :   missing - can we implement??
%       - zPosition     :   relative z position of microscope during data acquistion
%       - umPerPx_x     :   um per pixel conversion factor along x axis
%       - umPerPx_y     :   um per pixel conversion factor along y axis
%       - nCh           :   number of channels acquired
%       - channels      :   list of channels that are recorded (e.g. [2, 3])
%       - channelNames  :   list of corresponding channel names e.g. {Ch2, Ch3}
%       - channelColor  :   list of corresponding color for each channel e.g. {green, red}
%       - nBlocks       :   number of blocks
%       - nFrames       :   array of nFrames per block
%       - times         :   array of time vectors per recording/per block
%
%   see also getSciScanMetaData, getPrairieMetaData
%
% WB Eivind Hennestad, modified by Andreas Lande
% * Modified: AL removed blockfolders functionality and the redundant 
%   Prariescope support. AL also merged this function and the previous 
%   function called getSciScanMetaData() which was called from within this function.


%-- Open metadata file, Write in manually : laserPower, waveLength, fovCoordinates;
parts = strsplit(List.folder,'\'); % Clear"/roisignals" from path.
FolderMeta = strcat(parts{1},'\',parts{2},'\',parts{3},'\',parts{4},'\metadata');
List2 = dir(fullfile(FolderMeta,'*.mat'));
if ~isempty(List2)
    load(fullfile(FolderMeta,List2.name));
    metadata = meta2P;    
else
    %-- If metadata file is not in the motion corrected folder
    msgbox('Open raw recording metadata (configuration setting file), Write in manually : laserPower, waveLength, fovCoordinates');
    [fileName,filePath,~] = uigetfile('*.*');
    inifilepath = fullfile(filePath, fileName);
    inistring = fileread(inifilepath);
    metadata = struct();
    
    % Detect microscope used
    metadata.microscope = readVarIni(inistring,'microscope');
    if isempty(metadata.microscope) % Old recordings do not have this info	
        % Use the difference in root folder as indicator
        root_folder = readVarIni(inistring,'root.path');
        if root_folder(1) == 'E'
            metadata.microscope = 'OS2';
        elseif root_folder(1) == 'D'
            metadata.microscope = 'OS1';
        else
            metadata.microscope = 'OS1';
            warning('There is missing root path data. Could not find which microscope was used');
        end
    end

    % Get recording metadata
    metadata.xpixels=readVarIni(inistring,'x.pixels');
    metadata.ypixels=readVarIni(inistring,'y.pixels');
    aqusition_freq = readVarIni(inistring,'frames.p.sec');
    metadata.fps = aqusition_freq(1);
    metadata.dt = 1/aqusition_freq(1);

    % Detect settings for when the piezo is used.
    piz = readVarIni(inistring,'piezo.active');
    if (piz(2:5) == 'TRUE')
        metadata.piezoActive = true;
        metadata.piezoNumberOfPlanes = readVarIni(inistring,'frames.per.z.cycle');
        metadata.piezoImagingRateHz = readVarIni(inistring,'volume.rate.(in.Hz)');
        metadata.zDepthYm = (readVarIni(inistring,'z.spacing')*readVarIni(inistring,'no.of.planes'))-readVarIni(inistring,'z.spacing');

        % Detect if zig-zag or sawtooth mode for piezo is used
        piz_mode = readVarIni(inistring,'piezo.mode');
        if (piz_mode(2:5) == 'TRUE')
            metadata.piezoMode = 'saw';
        else
            metadata.piezoMode = 'zig';
        end

    else
        metadata.piezoActive = false;
    end

    % Add metadata.times. Should be a cell array of array with timevec.
    metadata.zoomFactor = readVarIni(inistring,'ZOOM');
    metadata.zPosition = abs(readVarIni(inistring,'setZ'));
    xfov = abs(readVarIni(inistring,'x.fov'));
    yfov = abs(readVarIni(inistring,'y.fov'));

    metadata.umPerPx_x = xfov / metadata.xpixels;
    metadata.umPerPx_y = yfov / metadata.ypixels;

    PMT1gain = readVarIni(inistring,'pmt1.gain');
    PMT2gain = readVarIni(inistring,'pmt2.gain');
    metadata.pmtGain = [PMT1gain, PMT2gain];

    %-- Get info about recorded channels
    % The PMT on OS1 and OS2 is switched
    if metadata.microscope == 'OS1'
        colors = {'Red', 'Green', 'N/A', 'N/A'};
    else
        colors = {'Green', 'Red', 'N/A', 'N/A'};
    end

    metadata.channelNumbers = [];
    metadata.channelNames = {};
    metadata.channelColor = {};

    for ch = 1:4
        if strcmp(strtrim(readVarIni(inistring, ['save.ch.', num2str(ch)])), 'TRUE')
            metadata.channelNames{end+1} = ['Ch', num2str(ch)];
            metadata.channelColor{end+1} = colors{ch};
            metadata.channelNumbers(end+1) = ch;
        end
    end

    metadata.nChannels = length(metadata.channelNumbers);

    %-- Get number of frames and estimate time onset for each
    metadata.nFrames = [];
    metadata.times = {};

    try
        metadata.nFrames(end+1) = readVarIni(inistring, 'no.of.frames.acquired');
    catch
        frame_count = readVarIni(inistring, 'frame.count');
        metadata.nFrames(end+1) = frame_count(1);
    end

    metadata.times{end+1} = linspace(0, metadata.nFrames(end)*metadata.dt, metadata.nFrames(end));
end

% Add fields for data that has to be manually added.
metadata.laserPower = [];
metadata.waveLength = [];
metadata.fovCoordinates = [];


end