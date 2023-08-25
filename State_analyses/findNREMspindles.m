function sData = findNREMspindles(sData)
% first sleep score session using episode marker function and detect
% spindles before running this function. 


% find NREM episodes > 30s
NREMep = sData.episodes.state == 'NREM' & sData.episodes.state_duration > 30;

% copy array in 2nd column
logidx = repmat(NREMep, 1, 2);

% convert table to array
NREMarr = table2array(sData.episodes(:,2:3));

% select NREM bouts 
NREM = NREMarr(logidx);

% reshape NREM start/end vector into original matrix and multiply by sample
% rate (2500)
NREMstartend = reshape(NREM, length(NREM)/2, 2).*2500;
sData.NREMepisodes = NREMstartend;

% initialize NREM spindle matrix
NREMspindles2 = false(size(sData.ephysdata2.absSpindleIdx));

% loop over nr of spindles
for i = 1:length(sData.ephysdata2.absSpindleIdx)
    
    % go through NREM bouts and check if spindle occurs in one of them
    aa = 0;
    t = 1;
    while aa < 1 && t <= length(NREMstartend)
        if sData.ephysdata2.absSpindleIdx(i) > NREMstartend(t,1) && sData.ephysdata2.absSpindleIdx(i) < NREMstartend(t,2)
           NREMspindles2(i,1) = true;
           aa = 1;
        end
        t = t +1;
    end
end

sData.ephysdata2.NREMAbsSpindleIdx = sData.ephysdata2.absSpindleIdx(NREMspindles2);

d1 = repmat(NREMspindles2,1,2);
d2 = sData.ephysdata2.spindleStartEnd(d1);
sData.ephysdata2.NREMspindleStartEnd = reshape(d2, length(d2)/2, 2);
sData.ephysdata2.NREMspindleSnips = sData.ephysdata2.spindleSnips(d1(:,1));

sData.ephysdata2.NREMspindleCycl = sData.ephysdata2.spindleCycl(NREMspindles2);
sData.ephysdata2.NREMdur = sData.ephysdata2.spindleDur(NREMspindles2);
sData.ephysdata2.NREMfreq = sData.ephysdata2.spindleFreq(NREMspindles2);
sData.ephysdata2.NREMamp = sData.ephysdata2.spindleAmp(NREMspindles2);

%% repeat for channel 3
NREMspindles3 = false(size(sData.ephysdata3.absSpindleIdx));

% loop over nr of spindles
for i = 1:length(sData.ephysdata3.absSpindleIdx)
    
    % loop over nr of NREM bouts
    aa = 0;
    t = 1;
    while aa < 1 && t <= length(NREMstartend)
        % check if spindle central time point occurs within a NREM bout 
        if sData.ephysdata3.absSpindleIdx(i) > NREMstartend(t,1) && sData.ephysdata3.absSpindleIdx(i) < NREMstartend(t,2)
           NREMspindles3(i,1) = true;
           aa = 1;
        end
        t = t +1;
    end
end

sData.ephysdata3.NREMAbsSpindleIdx = sData.ephysdata3.absSpindleIdx(NREMspindles3);

c1 = repmat(NREMspindles3,1,2);
c2 = sData.ephysdata3.spindleStartEnd(c1);
sData.ephysdata3.NREMspindleStartEnd = reshape(c2, length(c2)/2, 2);
sData.ephysdata3.NREMspindleSnips = sData.ephysdata3.spindleSnips(c1(:,1));

sData.ephysdata3.NREMspindleCycl = sData.ephysdata3.spindleCycl(NREMspindles3);
sData.ephysdata3.NREMdur = sData.ephysdata3.spindleDur(NREMspindles3);
sData.ephysdata3.NREMfreq = sData.ephysdata3.spindleFreq(NREMspindles3);
sData.ephysdata3.NREMamp = sData.ephysdata3.spindleAmp(NREMspindles3);

%% repeat for channel 4
NREMspindles4 = false(size(sData.ephysdata4.absSpindleIdx));

% loop over nr of spindles
for i = 1:length(sData.ephysdata4.absSpindleIdx)
    
    % loop over nr of NREM bouts
    aa = 0;
    t = 1;
    while aa < 1 && t <= length(NREMstartend)
        % check if spindle central time point occurs within a NREM bout 
        if sData.ephysdata4.absSpindleIdx(i) > NREMstartend(t,1) && sData.ephysdata4.absSpindleIdx(i) < NREMstartend(t,2)
           NREMspindles4(i,1) = true;
           aa = 1;
        end
        t = t +1;
    end
end

sData.ephysdata4.NREMAbsSpindleIdx = sData.ephysdata4.absSpindleIdx(NREMspindles4);

b1 = repmat(NREMspindles4,1,2);
b2 = sData.ephysdata4.spindleStartEnd(b1);
sData.ephysdata4.NREMspindleStartEnd = reshape(b2, length(b2)/2, 2);
sData.ephysdata4.NREMspindleSnips = sData.ephysdata4.spindleSnips(b1(:,1));

sData.ephysdata4.NREMspindleCycl = sData.ephysdata4.spindleCycl(NREMspindles4);
sData.ephysdata4.NREMdur = sData.ephysdata4.spindleDur(NREMspindles4);
sData.ephysdata4.NREMfreq = sData.ephysdata4.spindleFreq(NREMspindles4);
sData.ephysdata4.NREMamp = sData.ephysdata4.spindleAmp(NREMspindles4);

% calculate total NREM sleep minutes
totalNREM = sData.ephysdata2.lfp(sData.NREMepisodes(1,1):sData.NREMepisodes(1,2))';
for i = 2:length(sData.NREMepisodes)
    totalNREM = horzcat(totalNREM, sData.ephysdata2.lfp( sData.NREMepisodes(i,1):sData.NREMepisodes(i,2))');
end
sData.totalNREMmin = length(totalNREM)./ (2500*60);
% calculate spindle density (nr of spindles per NREM sleep minute)
sData.ephysdata2.NREMspindleDen = length(sData.ephysdata2.NREMAbsSpindleIdx)/sData.totalNREMmin;
sData.ephysdata3.NREMspindleDen = length(sData.ephysdata3.NREMAbsSpindleIdx)/sData.totalNREMmin;
sData.ephysdata4.NREMspindleDen = length(sData.ephysdata4.NREMAbsSpindleIdx)/sData.totalNREMmin;


