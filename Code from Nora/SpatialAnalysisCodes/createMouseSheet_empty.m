mouseID = 8007;

SavePath = 'C:\MATLAB\MOUSEINFO';
mouseInfo = struct();                     % Struct containing all information about the mouse

% MOUSEINFO
sData.mouseInfo.name = 'X';                                  % char      (REQUIRED) Mouse number according to lab convention
sData.mouseInfo.dateOfBirth = 'X.X.X';                     % char      (REQUIRED) yyyy.mm.dd
sData.mouseInfo.strain = 'X';               % char      (REQUIRED) Strain of animal used. Be specific. 
sData.mouseInfo.transgene = 'X';
sData.mouseInfo.Ordered = 'X';
sData.mouseInfo.sex = 'X';                                 % char      (REQUIRED) Male or female
sData.mouseInfo.SLCageNumber = 0;                           % double    (Optional) Science linker cage number
sData.mouseInfo.SLMouseNumber = 0;                          % double    (Optional) Science Linker animal number
sData.mouseInfo.SLMouseLitter = 0;
sData.mouseInfo.lightCycle = 'Light on 10pm-10am';              % char      (Optional) Light is on during the following time, ex: '10 am to 10 pm'.

sData.mouseInfo.surgeryDate = 'X.X.X';                          % char      (REQUIRED) yyyy.mm.dd
sData.mouseInfo.windowCoordinates = 'X';      % double	'AP -2.2 mm, ML 0 mm'    (REQUIRED?) [x, y] (mm) distance from bregma to center of window.
sData.mouseInfo.surgeryDoneBy = 'Nora';
sData.mouseInfo.surgeryProtocol = 'X';                 %'Virus injection + window implantation'  char      (Optional) should refer to a documented protocol
sData.mouseInfo.windowType = 'X';                 % char      Halfmoon/fullmoon
sData.mouseInfo.injectedVirusN1 = 'X';       % char      AAV1-CAG-flex-GCamp6s
sData.mouseInfo.injectedVirusN1Location = 'X';            %  'Left RSC'
sData.mouseInfo.injectedVirusN1NanoLPerSite = 0;                % double
sData.mouseInfo.injectedVirusN1NumberOfSites = 0;                % double
sData.mouseInfo.injectedVirusN1Depth = 0;                      % double (micron)
sData.mouseInfo.RecordedHemisphere = 'X';                     % char      (Optional) relevant for injections/recordings
%sData.mouseInfo.injectionCoordinates = '';                      % double    (Optional) nInj x 2 (um) from bregma

% save 
save(fullfile(SavePath,strcat('mouseInfo-',mouseID,'.mat')),'mouseInfo');
