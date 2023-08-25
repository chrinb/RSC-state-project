function singlets_cluster_data = swr_aligned_ecog_cluster(sData)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots mean WR-aligned ECoG activity for true SWR singlets 
% and different SWR clusters. 


%% Select SWRs for analysis
prompt = sprintf('All SWRs (1) | Rem. locom. SWRs (2) | Rem. clust. SWRs (3) | Rem- locom. & clust. SWRs (4) ');
swr_for_analysis = input(prompt);

RippleIdx     = sData.ephysdata.absRipIdx;

if swr_for_analysis == 2
    [RippleIdx,~] = riprun2(sData, RippleIdx);
elseif swr_for_analysis == 3
    [RippleIdx]   = RemoveRip(RippleIdx,[],[]);
elseif swr_for_analysis == 4
    [RippleIdx,~] = riprun2(sData, RippleIdx);
    [RippleIdx]   = RemoveRip(RippleIdx,[],[]);
end


% Load cluster data
swr_clust = swr_cluster_an(sData, RippleIdx);

% Load singlet data
[true_sing_idx, true_sing_nr] = true_singlet(sData, RippleIdx, []);

nr_of_seconds = 1;
RippleIdx = sData.ephysdata.absRipIdx(true_sing_idx);

% Create variables for analysis
n_ripples        = length(RippleIdx);
RSC_ECoG         = sData.ephysdata2.lfp;
win_length       = (nr_of_seconds*2500*2)+1;
time             = (-(nr_of_seconds*2500):(nr_of_seconds*2500))./2500;
sessionID        = sData.sessionInfo.sessionID;

% Specify whether to do baseline subtraction
prompt  = sprintf('Do baseline subtraction? (1) | everything else = no) ');
baseSub = input(prompt);
if baseSub ~= 1
    baseSub = [];
end

% Preallocate
singlet_ECoG_signals = zeros(true_sing_nr, length(time));

%% Run main analysis 

% Loop over true singlets
for nr_singlet = 1:true_sing_nr
    
    % find beginning and end of ECoG snippet centered on SWR
    singlet_signal_start = RippleIdx(nr_singlet) - (nr_of_seconds * 2500);
    singlet_signal_end   = RippleIdx(nr_singlet) + (nr_of_seconds * 2500);

    % skip SWRs occurring too early or late (i.e. the SWR-window
    % extends before or beyond recording length) 
    if singlet_signal_start > 1 && singlet_signal_end < length(RSC_ECoG)
        singlet_ECoG_signals(nr_singlet,:) = RSC_ECoG(singlet_signal_start:singlet_signal_end);
    end

    % Baseline subtraction
    if baseSub == 1
        % Baseline = mean activity in -1 to -0.8 sec before ripple peak
        baseline_vm                            = nanmean( singlet_ECoG_signals( nr_singlet,1:200));
        singlet_ECoG_signals(nr_singlet,:) = singlet_ECoG_signals( nr_singlet,:)-baseline_vm;
    end
    singlets_cluster_data{1,1} = 'True singlets';
    singlets_cluster_data{2,1} = singlet_ECoG_signals;
end

% Loop over SWR cluster types
for n_clust = 1:size(swr_clust,2)
    
    % Loop over nr of first SWRs in current cluster type
    for cluster_nr = 1:length(swr_clust{8, n_clust}  )
    
        cluster_signal_start = swr_clust{8, n_clust}(cluster_nr)  - (nr_of_seconds * 2500);
        cluster_signal_end   = swr_clust{8, n_clust}(cluster_nr)  + (nr_of_seconds * 2500);

        if cluster_signal_start > 1 && cluster_signal_end < length(RSC_ECoG)
            cluster_ECoG_signals(cluster_nr,:) = RSC_ECoG(cluster_signal_start:cluster_signal_end);
        end

        % Baseline subtraction
        if baseSub == 1
            % Baseline = mean activity in -1 to -0.8 sec before ripple peak
            baseline_vm                        = nanmean( cluster_ECoG_signals( cluster_nr,1:200));
            cluster_ECoG_signals(cluster_nr,:) = cluster_ECoG_signals( cluster_nr,:)-baseline_vm;
        end
    end
    singlets_cluster_data{1, n_clust+1} = n_clust+1;
    singlets_cluster_data{2, n_clust+1} = cluster_ECoG_signals;
    cluster_ECoG_signals = [];
end

figure, 

subplot(3,3, [1 4])
imagesc(singlets_cluster_data{2,1})


subplot(3,3, [2 5])
imagesc(singlets_cluster_data{2,2})

subplot(3,3, [3,6])
imagesc(singlets_cluster_data{2,3})

