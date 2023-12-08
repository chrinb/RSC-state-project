function sData = hierClust_Axons3(varargin)

% Written by ?. Modified by Christoffer Berge

% Function that merges correlated (axon) ROIs separately for multi-plane recordings.
% Identical function as "hierClust_Axons2"


%% Prepare different variables

% Load sData
sData                 = varargin{1,1};

% Get axon ROI indices
[~, ~, cell_idx_axon, ~, ~, all_idx_axon] = remove_cells_piezo(sData);

% Load DF/F
dff                   = sData.imdata.roiSignals(2).newdff;
smoothSignal          = 1;   
merged_ROI           = [];

% Determine how many planes contain axon ROIs
plane_idx = [1, 2, 3, 4];
plane_nr  = plane_idx( ~cellfun(@isempty, (cell_idx_axon)));

% Get 2P frame rate
imaging_sampling_rate = find_imaging_framerate(sData);

% Divide by 4 to get FPS per plane
imaging_sampling_rate = imaging_sampling_rate/4;

%% Smooth, correlate and cluster ROIs

sData.imdata.plane_indices_merged = [];
dff_merged                        = [];
dec_merged                        = [];
dff_merged_okada                  = [];

% Loop over planes
for current_plane = plane_nr
    
    temp_idx = sData.imdata.plane_indices(all_idx_axon) == current_plane;
    temp_dff = dff( all_idx_axon( temp_idx),:);

    % Remove NaNs from data
    nan_idx = isnan(temp_dff);
    temp_dff(nan_idx(:,1),:) = [];
    
    try 
        roiClustIDs = varargin{1,5};
    catch
        roiClustIDs = [];
    end
    
    % reduce noise effects
    if smoothSignal == 1
        dff_smooth = smoothdata(temp_dff,2,'gaussian',20);  
    end
    
    % Z-score data
    dff_zscore = zscore(dff_smooth, [],2);
    
    % Find nr of ROIs
    nRois       = size(dff_zscore,1);
    
    % Compute ROI signal correlations
    Y          = pdist(dff_zscore, 'correlation');
    roiVals    = squareform(Y);
    temp       = tril(roiVals, -1);
    % uniqueCorr = temp(temp ~= 0);
        
    % Cluster ROIs with correlation >= 0.7
    linkg  = 'average'; %settings.linkage.Selection;
    tree   = linkage(double(Y), linkg); %, 'correlation'); 
    T      = cluster(tree,'Cutoff', 0.4, 'Criterion', 'distance');
    nclust = max(T);
    
    % Plot cluster
    % figure()
    % H= dendrogram(tree, size(tree,1));
    % set(H,'LineWidth',1.5, 'Color', [0 0 0])
    % yline(0.4, 'LineWidth', 2, 'LineStyle','--')

    
    % sort the clusters by highest mean correlation
     nclustReal = []; 
     j          = 1;
    
    % if existing clustering does not exist, begin merging from scratch
    if isempty(roiClustIDs)
        
        % loop over nr of clusters
        for nr_clust=1:nclust   
        
            % for all clusters containing more than 1 ROI
            if length(find(T== nr_clust)) > 1
        
                % Average the values in the correlation distance matrix for those
                % ROIs
                nreal(j)   = mean( mean( roiVals(T==nr_clust, T==nr_clust),'omitnan' ), 'omitnan');
                nclustReal = [nclustReal nr_clust];
                j = j+1;
            end
        end
        
        % check if some ROIs have been clustered. If not, skip remaining
        % analysis
        if exist('nreal')
            [qsort_clust, indclust] = sort(nreal,'descend');
            
            finalClust = indclust(qsort_clust < 0.3);
            finalClust = nclustReal(finalClust);
            
            mean(qsort_clust ~= 1)
            
            % merge axons within highest correlated clusters
            T(:,2)    = 1:nRois;
            T         = sortrows(T,1);
            clustRois = [];
        
            % Check if there are clusters
            if isempty(finalClust)
                disp('No significantly correlated clusters!')
                   merged_ROI        = temp_dff;
                   roiClustIDs = [];
            % If yes, loop over nr of clusters    
            else
                t = 1;
                for i = 1:length(finalClust)
                    
                    % For each cluster, find the indicies of clustered ROIs  
                    A = find(T(:,1) == finalClust(i));
                    % Check that there are more than 1 ROI in cluster
                    if length(A) > 1
                        rois                     = T(A,2);
                        clustRois                = [clustRois; rois];
                        roiClustIDs(i).rois      = rois; 
                        clustSignals(t,:)        = mean(temp_dff(rois, :), 'omitnan');
                        t = t + 1;        
                    end
                    A = [];
                    
                end 
                    temp_dff(clustRois,:)   = []; 
                    merged_ROI         = [temp_dff;        clustSignals];    
                    n_unchanged_rois   = size(temp_dff,1);
                    cluster_new_roi_id = [size(temp_dff,1)+1:( size(temp_dff,1)+size(clustSignals,1)); 1:size(clustSignals,1)];
            end
    
            else
                   disp('No significantly correlated clusters!')
                   merged_ROI        = temp_dff;
    %                merged_ROI_deconv = deconv;
                   roiClustIDs = [];
        end
    
    % If existing roi clustering IDs is given as input argument, apply this to 
    % current sData
    else
       
        % Loop over nr of ROI clusters
        for i = 1:size(roiClustIDs.rois,2)
            
            % Get indicies of clustered ROIs in each cluster
            rois = roiClustIDs(i).rois;
            
            
            % Average clustered ROIs
            clustSignals        = mean(temp_dff(rois,:), 'omitnan');
    
            % Remove clustered ROIs
            temp_dff(rois,:)        = [];
            temp_dff_zscore(rois,:) = [];
    %         deconv(rois,:)     = [];
        
            merged_ROI        = [temp_dff;        clustSignals];
        end
    end
    
    % Filter merged axonal DF/F signals using Okada filter
    merged_ROI_okada = okada(merged_ROI,2);
    
    % Deconvolve merged DF/F
    P = nansen.twophoton.roisignals.getDeconvolutionParameters();
    P.modelType  = 'autoar';
    P.spikeSnr   = 2;
    P.sampleRate = imaging_sampling_rate;
    [dec, ~, ~]  = nansen.twophoton.roisignals.deconvolveDff(merged_ROI, P);


    % Concatenate merged ROIs from different planes into one array
    dff_merged                        = [dff_merged; merged_ROI];
    dff_merged_okada                  = [dff_merged_okada; merged_ROI_okada];
    dec_merged                        = [dec_merged; dec];   
    sData.imdata.plane_indices_merged = [sData.imdata.plane_indices_merged; ones(size(merged_ROI,1),1) * current_plane];
    merged_ROI                        = [];
    merged_ROI_okada                  = [];
    dec                               = [];
    nreal                             = [];
end
%% Store output in sData
sData.imdata.roiSignals(2).mergedAxonsDff     = dff_merged;
sData.imdata.roiSignals(2).mergedAxonsDffFilt = dff_merged_okada;
sData.imdata.roiSignals(2).mergedAxonsDec     = dec_merged;
sData.analysis.roiClustIDs                    = roiClustIDs; % Save IDs of clustered ROIs

try
    sData.analysis.n_unchanged_rois               = n_unchanged_rois; %
    sData.analysis.clusterID                      = cluster_new_roi_id;
catch
end
