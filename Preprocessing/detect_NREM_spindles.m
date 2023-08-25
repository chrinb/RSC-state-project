function sData = detect_NREM_spindles(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that finds sleep spindles occurring inside NREM bouts > 30s.
% User first have to sleep score the session using "episode marker" function and detect
% spindles before running this function. NOTE: same function as "findNREMspindles" 
% except that this only analyzes one channel

sData = varargin{1,1};

if length(varargin) < 2
    prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
    spindle_band_select = input(prompt);

    if spindle_band_select == 1
        spindle_select = [];
    elseif spindle_band_select == 2
        spindle_select = '1016';
    end
else
    spindle_select = varargin{1,2};
end

% find NREM episodes > 30s (if there are no "episodeMerge" field, use the
% "episode" field in sData instead).
try
    NREMep       = sData.episodeMerge.state == 'NREM' & sData.episodeMerge.state_duration > 30;
    NREMarr      = table2array(sData.episodeMerge(NREMep,2:3));
catch
    NREMep       = sData.episodes.state == 'NREM' & sData.episodes.state_duration > 30;
    NREMarr      = table2array(sData.episodes(NREMep,2:3));
end

NREMstartend = NREMarr*2500;
sData.NREMepisodes = NREMstartend;             

% initialize NREM spindle matrix
str_spin_idx = strcat('absSpindleIdx', spindle_select);
NREMspindles2 = false(size(sData.ephysdata2.(str_spin_idx)) );

% loop over nr of spindles
for i = 1:length(sData.ephysdata2.(str_spin_idx))
    
    % go through NREM bouts and check if spindle occurs in one of them
    aa = 0;
    t = 1;
    while aa < 1 && t <= size(NREMstartend,1)
        if sData.ephysdata2.(str_spin_idx)(i) > NREMstartend(t,1) && sData.ephysdata2.(str_spin_idx)(i) < NREMstartend(t,2)
           NREMspindles2(i,1) = true;
           aa = 1;
        end
        t = t +1;
    end
end

NREM_spin_idx_str = strcat('NREMAbsSpindleIdx', spindle_select);
NREM_start_end_str = strcat('NREMspindleStartEnd', spindle_select);
NREM_spin_snip_str = strcat('NREMspindleSnips', spindle_select);
NREM_spin_cycl_str = strcat('NREMspindleCycl', spindle_select);
NREM_spin_dur_str = strcat('NREMdur', spindle_select);
NREM_spin_freq_str = strcat('NREMfreq', spindle_select);
NREM_spin_amp_str = strcat('NREMamp', spindle_select);

spin_start_end_str = strcat('spindleStartEnd', spindle_select);
spin_snip_str = strcat('spindleSnips', spindle_select);
spin_cycl_str = strcat('spindleCycl', spindle_select);
spin_dur_str = strcat('spindleDur', spindle_select);
spin_freq_str = strcat('spindleFreq', spindle_select);
spin_amp_str = strcat('spindleAmp', spindle_select);

sData.ephysdata2.(NREM_spin_idx_str) = sData.ephysdata2.(str_spin_idx)(NREMspindles2);

d1 = repmat(NREMspindles2,1,2);
d2 = sData.ephysdata2.(spin_start_end_str)(d1);
sData.ephysdata2.(NREM_start_end_str) = reshape(d2, length(d2)/2, 2);
sData.ephysdata2.(NREM_spin_snip_str) = sData.ephysdata2.(spin_snip_str)(d1(:,1));

sData.ephysdata2.(NREM_spin_cycl_str) = sData.ephysdata2.(spin_cycl_str)(NREMspindles2);
sData.ephysdata2.(NREM_spin_dur_str) = sData.ephysdata2.(spin_dur_str)(NREMspindles2);
sData.ephysdata2.(NREM_spin_freq_str) = sData.ephysdata2.(spin_freq_str)(NREMspindles2);
sData.ephysdata2.(NREM_spin_amp_str) = sData.ephysdata2.(spin_amp_str)(NREMspindles2);

%%  Calculate total NREM sleep minutes

% Extract first NREM episode
totalNREM = sData.ephysdata2.lfp(sData.NREMepisodes(1,1):sData.NREMepisodes(1,2))';

% Check if there are additional NREM episode in session. If yes, extract
% and concatenate episodes.
if size(sData.NREMepisodes,1) > 1

    % Loop over NREM episodes starting from the second (first is already
    % extracted)
    for i = 2:length(sData.NREMepisodes)
        % Concatenate episodes
        totalNREM = horzcat(totalNREM, sData.ephysdata2.lfp( sData.NREMepisodes(i,1):sData.NREMepisodes(i,2))');
    end
end

% calculate total REM sleep minutes
% totalREM = sData.ephysdata2.lfp(sData.NREMepisodes(1,1):sData.NREMepisodes(1,2))';
% for i = 2:length(sData.NREMepisodes)
%     totalNREM = horzcat(totalNREM, sData.ephysdata2.lfp( sData.NREMepisodes(i,1):sData.NREMepisodes(i,2))');
% end

sData.totalNREMmin = length(totalNREM)./ (2500*60);
% sData.totalREMmin  = length(totalREM)

% calculate spindle density (nr of spindles per NREM sleep minute)
sData.ephysdata2.NREMspindleDen = length(sData.ephysdata2.NREMAbsSpindleIdx)/sData.totalNREMmin;


