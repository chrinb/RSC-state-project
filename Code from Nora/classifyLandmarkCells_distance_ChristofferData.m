

file  ={'m6152-20230201-021','m6152-20230202-025','m6152-20230203-028',...
   'm6157-20230515-023','m6159-20230515-022','m6161-20230513-015','m6163-20230515-021'};

% position of cues:
% cue1 =  [42 47 122 128]./2;
% cue2  = [63 68 100 106]./2;


% cue 1 distance : 60 cm /.2 cm
% cue 2 distance : 20 cm /.2 cm
binGapsPos = 1:1.5:157;
cueIdx = [1 1 2 2];

for i = 1:length(file)
    loaddrct =fullfile('/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Andreas/Christopher/CueTuned/',file{i});
    savedrct =fullfile('/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Andreas/Christopher/CueCells_pm1/',file{i});
    if ~exist(savedrct); mkdir(savedrct); end
    
    load(fullfile('/Users/annachristinagarvert/UIO Physiology Dropbox Dropbox/Lab Data/Ann Christin Garvert/Andreas/Christopher/',file{i},'/sData.mat'))
    
    load(fullfile(loaddrct,'rmaps.mat'))
    load(fullfile(loaddrct,'rmaps_dff.mat'))
        
    
    params.first_landmark_bins = floor(sData.imdata.cues.C1A/1.5)-1:ceil(sData.imdata.cues.C1B/1.5)+1;
    params.second_landmark_bins = floor(sData.imdata.cues.C4A/1.5)-1:ceil(sData.imdata.cues.C4B/1.5)+1;

    lcsA= sb.classify.landmarkCellsByShuffling(rmaps,params);
    
    params.first_landmark_bins = floor(sData.imdata.cues.C2A/1.5)-1:ceil(sData.imdata.cues.C2B/1.5)+1;
    params.second_landmark_bins = floor(sData.imdata.cues.C3A/1.5)-1:ceil(sData.imdata.cues.C3B/1.5)+1;

    lcsB= sb.classify.landmarkCellsByShuffling(rmaps,params);

    save(fullfile(savedrct,'lcsA.mat'),'lcsA')
    save(fullfile(savedrct,'lcsB.mat'),'lcsB')


    mkdir(fullfile(savedrct,'lcsA'))
    mkdir(fullfile(savedrct,'lcsB'))


        
end
    




