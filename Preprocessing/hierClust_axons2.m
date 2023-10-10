function sData = hierClust_axons2(varargin)

% Written by ?. Modified by Christoffer Berge

% Merge ROIs that have highly correlated signals, e.g. segments of the 
% same axon or dendrite using multiple, concatenated sessions as input.

% THIS IS THE CODE TO USE NOW! (08.12.22)

%% Get ROI indices

if nargin == 1 
    sData        = varargin{1,1};
    [pc_rois, ~] = remove_cells(sData);
    dff          = sData.imdata.roiSignals(2).newdff(pc_rois,:);
    smoothSignal = 1;   
%     decmonv        = sData.imdata.roiSignals(2).ciaDeconvolved(pc_rois,:);
else
    dff          = varargin{1,1};
    smoothSignal = varargin{1,2};
%     deconv       = varargin{1,3};
end

merged_ROI        = [];

%% Smooth, correlate and cluster ROIs

% Remove NaNs from data
nan_idx = isnan(dff);
dff(nan_idx(:,1),:) = [];

try 
    roiClustIDs = varargin{1,5};
catch
    roiClustIDs = [];
end

% reduce noise effects
if smoothSignal == 1
    dff_smooth = smoothdata(dff,2,'gaussian',20);  
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
T      = cluster(tree,'Cutoff', 0.3, 'Criterion', 'distance');
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
               merged_ROI        = dff;
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
                    clustSignals(t,:)        = mean(dff(rois, :), 'omitnan');
                    t = t + 1;        
                end
                A = [];
                
            end 
                dff(clustRois,:)   = []; 
                merged_ROI         = [dff;        clustSignals];    
                n_unchanged_rois   = size(dff,1);
                cluster_new_roi_id = [size(dff,1)+1:( size(dff,1)+size(clustSignals,1)); 1:size(clustSignals,1)];
        end

        else
               disp('No significantly correlated clusters!')
               merged_ROI        = dff;
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
        clustSignals        = mean(dff(rois,:), 'omitnan');

        % Remove clustered ROIs
        dff(rois,:)        = [];
        dff_zscore(rois,:) = [];
%         deconv(rois,:)     = [];
    
        merged_ROI        = [dff;        clustSignals];
    end
end

% Filter merged axonal DF/F signals using Okada filter
merged_ROI_okada = okada(merged_ROI,2);

% Deconvolve merged DF/F
P = nansen.twophoton.roisignals.getDeconvolutionParameters();
P.modelType = 'autoar';
P.spikeSnr  = 2;
[dec, ~, ~] = nansen.twophoton.roisignals.deconvolveDff(merged_ROI, P);
%% Store output in sData
sData.imdata.roiSignals(2).mergedAxonsDff     = merged_ROI;
sData.imdata.roiSignals(2).mergedAxonsDffFilt = merged_ROI_okada;
sData.imdata.roiSignals(2).mergedAxonsDec     = dec;
sData.imdata.roiClustIDs                    = roiClustIDs; % Save IDs of clustered ROIs
try
    sData.imdata.n_unchanged_rois               = n_unchanged_rois; %
    sData.imdata.clusterID                      = cluster_new_roi_id;
catch
end
