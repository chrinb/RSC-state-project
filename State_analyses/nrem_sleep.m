function NREM_episodes = nrem_sleep(varargin)

% Written by Christoffer Berge | Vervaeke Lab

sData = varargin{1,1};

% Function that finds NREM sleep bouts in recording from sData.episodes
% field after sleep scoring using Begonia. 

% If only input argument is sData, find all un-processed NREM episodes
if nargin == 1
    NREMep = sData.episodes.state == 'NREM';
    % check that there is NREM episodes in recording
    if ~isempty(NREMep)
        NREMarr       = table2array(sData.episodes(NREMep,2:3));
        NREM_episodes = NREMarr*2500;
    end

% If a second input argument is given, extract the merged NREM episodes
elseif nargin == 2

    NREMep = sData.episodeMerge.state == 'NREM';
    % check that there is NREM episodes in recording
    if ~isempty(NREMep)
        NREMarr       = table2array(sData.episodeMerge(NREMep,2:3));
        NREM_episodes = NREMarr*2500;
    end
end
