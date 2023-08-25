function output = plot_swr_modulation_nansen(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Plot SWR modulated cells 

% Select particular dataset
prompt = sprintf('Specify substring for particular dataset? ');
str   = input(prompt, 's');
if ~isempty(str)
    varargin = [varargin, str]; 

end

% Load data and mouse IDs (user types in e.g. "m6120")
[data, miceIDs] = load_multi_session_data(varargin);

% Find nr of mice in data
n_mice = numel(fieldnames(data));

mouse_data    = cell(n_mice,7);
mouse_data_SE = cell(n_mice,2);
% Loop over nr of mice
n_rois            = 0;     
n_activated_rois  = 0;
n_suppressed_rois = 0;
mean_swr_mod_act = [];
mean_swr_mod_sup = [];
roi_act_idx = [];
roi_sup_idx = [];
all_mice_n_rois_act = 0;
all_mice_n_rois_sup = 0;
all_mice_n_rois     = 0;
all_mice_mean_act   = [];
all_mice_mean_sup   = [];
all_mice_roi_act    = [];
all_mice_roi_sup    = [];
% Select which datasets to extract
if strcmp(varargin{1,1}, 'dff')
    idx1 = 2;
    idx2 = 7;
elseif strcmp(varargin{1,1}, 'deconv')
    idx1 = 4;
    idx2 = 10;
end

for mouse_nr = 1:n_mice

    % Get all data for current mouse
    temp_data = data.(miceIDs{mouse_nr}).mean_swr_mod;
%     temp_rois
    % loop over nr of sessions
    n_rois = 0;
    n_activated_rois = 0;
    n_suppressed_rois = 0;
    for session_n = 1:size(temp_data,1)

        % if only one session put data into cell
%         if size(temp_data,1) == 1
%             temp_cell = cell(1,1);
%             temp_cell{1,1} = temp_data;
%             temp_data = temp_cell;
%         end
    
        % get rois in current session
        temp_rois        = temp_data{session_n, 1}.rois;  
        temp_roi_act_idx = temp_data{session_n, 1}.sorted_rois{1, 5};
        temp_roi_sup_idx = temp_data{session_n, 1}.sorted_rois{1, 10};
        temp_data_act    = temp_data{session_n, 1}.sorted_rois{1, idx1};
        temp_data_sup    = temp_data{session_n, 1}.sorted_rois{1, idx2}; 

        % get overall nr of ROIs
        n_rois           = n_rois + numel(temp_rois);
        n_activated_rois  = n_activated_rois + numel(temp_data{session_n, 1}.sorted_rois{1, 5}    );
        n_suppressed_rois = n_suppressed_rois + numel( temp_data{session_n, 1}.sorted_rois{1, 10}   );

        % get ROI ids
        roi_act_idx = [roi_act_idx; temp_roi_act_idx];
        roi_sup_idx = [roi_sup_idx; temp_roi_sup_idx];
        % get mean data

        mean_swr_mod_act = [mean_swr_mod_act; temp_data_act];
        mean_swr_mod_sup = [mean_swr_mod_sup; temp_data_sup];

    end
    mouse_data{mouse_nr, 1}    = mean_swr_mod_act;
    mouse_data{mouse_nr, 2}    = mean_swr_mod_sup;
    mouse_data{mouse_nr, 3}    = n_rois;
    mouse_data{mouse_nr, 4}    = n_activated_rois;
    mouse_data{mouse_nr, 5}    = n_suppressed_rois;
    mouse_data{mouse_nr, 6}    = roi_act_idx;
    mouse_data{mouse_nr, 7}    = roi_sup_idx;

    mouse_data_SE{mouse_nr, 1} = std(mean_swr_mod_act, 'omitnan')./sqrt( size(mean_swr_mod_act, 1));
    mouse_data_SE{mouse_nr, 2} = std(mean_swr_mod_sup, 'omitnan')./sqrt( size(mean_swr_mod_sup, 1));

    
    all_mice_n_rois_act = all_mice_n_rois_act + mouse_data{mouse_nr, 4};
    all_mice_n_rois_sup = all_mice_n_rois_sup + mouse_data{mouse_nr, 5} ;
    all_mice_n_rois     = all_mice_n_rois + mouse_data{mouse_nr, 3};
    all_mice_mean_act   = [all_mice_mean_act; mouse_data{mouse_nr, 1}];
    all_mice_mean_sup   = [all_mice_mean_sup; mouse_data{mouse_nr, 2}];
    all_mice_roi_act    = [all_mice_roi_act; mouse_data{mouse_nr, 6}];
    all_mice_roi_sup    = [all_mice_roi_sup; mouse_data{mouse_nr, 7}];
    
end

all_mice_mean_actSE = std(all_mice_mean_act, 'omitnan')./sqrt( size(all_mice_mean_act, 1));
all_mice_mean_supSE = std(all_mice_mean_sup, 'omitnan')./sqrt( size(all_mice_mean_sup, 1));

% Compute percentage modulated cells
frac_activated     = all_mice_n_rois_act/all_mice_n_rois;
frac_suppressed    = all_mice_n_rois_sup/all_mice_n_rois;
frac_not_modulated = (all_mice_n_rois - (all_mice_n_rois_act+all_mice_n_rois_sup)) /all_mice_n_rois;

x = [ frac_activated, frac_suppressed, frac_not_modulated];

% Plot data
figure,
p = pie(x);

pText         = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt           = {'Activated ';'Suppressed ';'Not modulated '}; 
combinedtxt   = strcat(txt,percentValues); 

legend(txt,'Location','southoutside','Orientation','horizontal')

output = struct();
output.nrois      = all_mice_n_rois;
output.act_rois   = all_mice_n_rois_act;
output.mean_act   = all_mice_mean_act;
output.mean_actSE = all_mice_mean_actSE;
output.sup_rois   = all_mice_n_rois_sup;
output.mean_sup   = all_mice_mean_sup;
output_mean_supSE = all_mice_mean_supSE;
% ouput.
% HAVE TO FIX THIS
%     mouse_data{mouse_nr, 3}    = n_rois;
%     mouse_data{mouse_nr, 4}    = n_activated_rois;
%     mouse_data{mouse_nr, 5}    = n_suppressed_rois;
%     mouse_data{mouse_nr, 6}    = roi_act_idx;
%     mouse_data{mouse_nr, 7}    = roi_sup_idx;