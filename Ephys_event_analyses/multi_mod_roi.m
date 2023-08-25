function mod_roi = multi_mod_roi(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Code that allows user to get the mean DF/F during slow wave activity for multiple
% sessions and store data in a struct. For the first session to be analyzed
% user only inputs "sData". For every following session add the output (e.g.,
% "slow_mean_signal" as a second argument and the code will update the struct 
% with the data from that sessions without overwriting the previous data. 

sData = varargin{1,1};
if length(varargin) > 1
    mod_roi = varargin{1,2};
end

prompt = sprintf('SWR (1) or spindle (2)? ');
select = input(prompt);

%% Locate file name in folder (for saving file)
folder_name = dir;
folder_cell = strfind({folder_name.name}, '.mat');
for i = 1:length(folder_cell)
    if folder_cell{1,i} > 0
    idx = i;
    end
end

if length(idx) > 1
    error('Idx is > 1, could not locate sData file name')
end

fileName = folder_name(idx).name;

%% Run slow wave analysis

if select == 1
    mod_event =  ripple_mod(sData);
elseif select == 2
    mod_event = spindle_mod(sData);
end


% for sleep recordings
names = {'_ROIs_activated' , '_ROIs_activated_idx', '_ROIs_suppressed',...
        '_ROIs_suppressed_idx', '_activated_sig_to_plot', '_suppressed_sig_to_plot'...
        '_ROIs_unclassified', '_ROIs_unclassified_idx'};
fileName = matlab.lang.makeValidName(fileName);
fileNames = strcat(fileName, names);
for i = 1:length(names)
    mod_roi.(fileNames{i}) = mod_event{1,i};
end

clearvars -except mod_roi
