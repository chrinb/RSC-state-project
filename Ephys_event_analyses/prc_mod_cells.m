function combinedtxt = prc_mod_cells(rois_for_an, sorted_rois)


% Written by Christoffer Berge | Vervaeke Lab

% Function that calculates the percentage of modulated ROIs in session. 

total_n_rois    = size(rois_for_an,2);
activated_rois  = [];
try
    activated_rois = size(sorted_rois{1, 5},1 );
catch
end

suppressed_rois = [];
try 
    suppressed_rois = size(sorted_rois{1,10},1);
catch
end

% Compute percentage modulated cells
frac_activated     = activated_rois/total_n_rois;
frac_suppressed    = suppressed_rois/total_n_rois;
frac_not_modulated = (total_n_rois - (activated_rois+suppressed_rois)) /total_n_rois;

x = [ frac_activated, frac_suppressed, frac_not_modulated];

% Plot data
figure,
p = pie(x);

pText         = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt           = {'Activated ';'Suppressed ';'Not modulated '}; 
combinedtxt   = strcat(txt,percentValues); 

legend(txt,'Location','southoutside','Orientation','horizontal')

% pText(1).String = combinedtxt(1);
% pText(2).String = combinedtxt(2);
% pText(3).String = combinedtxt(3);


% Adjust label position
% p(4).Position = [1 -1 0];
% p(2).Position = [-.4 1.1 1];
% p(6).Position = [0.4 1.1 1];