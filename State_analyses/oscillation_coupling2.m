function output = oscillation_coupling2(varargin)

% Function that looks for sequential slow oscillation & spindles, and
% delta waves & spindles. Based on method outlined in Kim et al. (2019)
% "Competing Roles of Slow Oscillations and Delta Waves in Memory 
% Consolidation versus Forgetting". 

sData                        = varargin{1,1};
signal                       = sData.ephysdata2.lfp;
signal_filt                  = sData.ephysdata2.deltaband;
hpc_lfp                      = sData.ephysdata.lfp;
time                         = (0:length(signal)-1)/2500;
srate                        = 2500;
ECoG_spindle                 = sData.ephysdata2.NREMspindleStartEnd/2500;
spindle_center               = sData.ephysdata2.NREMAbsSpindleIdx;

[nrem_SO, nrem_delta_waves] = slow_wave_an(sData);

%% SO-Spindle nesting analysis
closest_so_idx = zeros( length(spindle_center),1);

for k = 1:length(spindle_center)
    spindle = spindle_center(k);
        [~, closest_so_idx(k)] = min( abs(spindle - nrem_SO(:,3)));
end

% Define nesting time window
nesting_win_end   = srate*1;
nesting_win_start = srate*.1;

% Find nested SO-spindles where the time difference between the peak of the
% spindle and up-state of the linked SO is between -0.5s and 1.0s.
putative_nested_SO = nrem_SO(closest_so_idx, 3);

nested_SO = [];
for k = 1:length(spindle_center)
    
    % only include SOs that occur before spindle
    if spindle_center(k) - putative_nested_SO(k) > 0
        if spindle_center(k) - putative_nested_SO(k) > nesting_win_start && ...
                spindle_center(k) - putative_nested_SO(k) < nesting_win_end
            nested_SO = vertcat( nested_SO, putative_nested_SO(k));
        end
    end
end


nested_spindle = [];
for k = 1:length(putative_nested_SO)
    
    % only include SOs that occur before spindle
    if  putative_nested_SO(k) - spindle_center(k) < 0
        if abs( putative_nested_SO(k) - spindle_center(k) ) > nesting_win_start && ...
               abs( putative_nested_SO(k) - spindle_center(k) ) < nesting_win_end
            nested_spindle = vertcat( nested_spindle, spindle_center(k));
        end
    end
end


spindle_idx        = ismember(spindle_center, nested_spindle );
SO_idx             = ismember(nrem_SO(:,3), nested_SO);
enumerate_so       = 1:size(nrem_SO,1);
enumerate_spindles = 1:length(spindle_center);

%% Delta-spindle nesting analysis
closest_delta_idx = zeros( length(spindle_center),1);

for k = 1:length(spindle_center)
    spindle = spindle_center(k);
        [~, closest_delta_idx(k)] = min( abs(spindle - nrem_delta_waves(:,3)));
end
putative_nested_delta = nrem_delta_waves(closest_delta_idx, 3);


nested_delta = [];
for k = 1:length(spindle_center)
    
    % only include SOs that occur before spindle
    if spindle_center(k) - putative_nested_delta(k) > 0
        if spindle_center(k) - putative_nested_delta(k) > nesting_win_start && ...
                spindle_center(k) - putative_nested_delta(k) < nesting_win_end
            nested_delta = vertcat( nested_delta, putative_nested_delta(k));
        end
    end
end


nested_spindle_delta = [];
for k = 1:length(putative_nested_delta)
    
    % only include SOs that occur before spindle
    if  putative_nested_delta(k) - spindle_center(k) < 0
        if abs( putative_nested_delta(k) - spindle_center(k) ) > nesting_win_start && ...
               abs( putative_nested_delta(k) - spindle_center(k) ) < nesting_win_end
            nested_spindle_delta = vertcat( nested_spindle_delta, spindle_center(k));
        end
    end
end


spindle_delta_idx        = ismember(spindle_center, nested_spindle_delta );
delta_idx             = ismember(nrem_delta_waves(:,3), nested_delta);
enumerate_delta       = 1:size(nrem_delta_waves,1);
%% Plot
figure,
plot(time,signal_filt)
hold on
plot(time,signal*.4-.2)

for i = enumerate_spindles(spindle_idx)
    z = [ECoG_spindle(i,1) ECoG_spindle(i,1) ECoG_spindle(i,2) ECoG_spindle(i,2) ];
    v = [-.5 .5 .5 -.5];
    patch(z, v, 'red', 'edgecolor', 'none', 'FaceAlpha', .3);
end

for i = enumerate_so(SO_idx)
    x = [nrem_SO(i,1) nrem_SO(i,1) nrem_SO(i,3) nrem_SO(i,3)]/2500;
    y = [-.5 .5 .5 -.5];
    patch(x, y, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .3);
end
set(gca,'ylim',[-.5 .5])


figure,
plot(time,signal_filt)
hold on
plot(time,signal*.4-.2)

for i = enumerate_spindles(spindle_delta_idx)
    z = [ECoG_spindle(i,1) ECoG_spindle(i,1) ECoG_spindle(i,2) ECoG_spindle(i,2) ];
    v = [-.5 .5 .5 -.5];
    patch(z, v, 'red', 'edgecolor', 'none', 'FaceAlpha', .3);
end

for i = enumerate_delta(delta_idx)
    x = [nrem_delta_waves(i,1) nrem_delta_waves(i,1) nrem_delta_waves(i,3) nrem_delta_waves(i,3)]/2500;
    y = [-.5 .5 .5 -.5];
    patch(x, y, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .3);
end
set(gca,'ylim',[-.5 .5])