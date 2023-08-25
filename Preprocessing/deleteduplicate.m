function sData = deleteduplicate(sData)

%Written by Christoffer Berge | Vervaeke Lab

prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
spindle_band_select = input(prompt);

if spindle_band_select == 1
    spindle_select = [];
elseif spindle_band_select == 2
    spindle_select = '1016';
end


% Delete duplicate spindles
spin_idx_str = strcat('absSpindleIdx', spindle_select);
spin_start_end_str = strcat('spindleStartEnd', spindle_select);
spin_snip_str = strcat('spindleSnips', spindle_select);
spin_cycl_str = strcat('spindleCycl', spindle_select);
spin_dur_str = strcat('spindleDur', spindle_select);
spin_freq_str = strcat('spindleFreq', spindle_select);
spin_amp_str = strcat('spindleAmp', spindle_select);

spindleLocs = sData.ephysdata2.(spin_idx_str);
spindleSnips = sData.ephysdata2.(spin_snip_str);
spindleStartEnd = sData.ephysdata2.(spin_start_end_str);
spindleCycl = sData.ephysdata2.(spin_cycl_str);
spindleDur = sData.ephysdata2.(spin_dur_str);
spindleFreq = sData.ephysdata2.(spin_freq_str);
spindleAmp = sData.ephysdata2.(spin_amp_str);

[~, doubleLocs] = unique(spindleLocs);
x = 1:length(spindleLocs);
x(doubleLocs) = [];

if length(x) > 1
    for i = 1:length(x)
        i = x(i);
        spindleLocs(i) = 0;
        spindleSnips(i).lfp = NaN(20001,1);
        spindleStartEnd(i,:) = 0;
        spindleCycl(i) = 0;
        spindleDur(i) = 0;
        spindleFreq(i) = 0;
        spindleAmp(i) = 0;
    end
    % delete duplicate spindles
    spindleLocs(spindleLocs == 0) =[];
    %spindleSnips(spindleSnips == 0) =[];
    spindleStartEnd( ~any(spindleStartEnd,2), : ) = [];
    spindleCycl(spindleCycl == 0) =[];
    spindleDur(spindleDur == 0) =[];
    spindleFreq(spindleFreq == 0) =[];
    spindleAmp(spindleAmp == 0) =[];
        
    sData.ephysdata2.(spin_idx_str)       = spindleLocs;
    sData.ephysdata2.(spin_snip_str)      = spindleSnips;  
    sData.ephysdata2.(spin_start_end_str) = spindleStartEnd;
    sData.ephysdata2.(spin_cycl_str)      = spindleCycl;
    sData.ephysdata2.(spin_dur_str)       = spindleDur;
    sData.ephysdata2.(spin_freq_str)      = spindleFreq;
    sData.ephysdata2.(spin_amp_str)       = spindleAmp;
elseif length(x) == 1
    spindleLocs(x) = [];
    spindleSnips(x).lfp = NaN(20001,1);
    spindleStartEnd(x,:) = [];
    spindleCycl(x) = [];
    spindleDur(x) = [];
    spindleFreq(x) = [];
    spindleAmp(x) = [];
    
    sData.ephysdata2.(spin_idx_str)       = spindleLocs;
    sData.ephysdata2.(spin_snip_str)      = spindleSnips;  
    sData.ephysdata2.(spin_start_end_str) = spindleStartEnd;
    sData.ephysdata2.(spin_cycl_str)      = spindleCycl;
    sData.ephysdata2.(spin_dur_str)       = spindleDur;
    sData.ephysdata2.(spin_freq_str)      = spindleFreq;
    sData.ephysdata2.(spin_amp_str)       = spindleAmp;
end
