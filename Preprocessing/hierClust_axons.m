function [merged_ROI, merged_ROI_deconv, roiClustIDs] = hierClust_axons(varargin)

%PURPOSE: merge ROIs that have highly correlated signals, e.g. segments of
%         the same axon or dendrite

%INPUT: ROImtx, an n x m matrix where n is number of deconvolved ROIs, and m is
%       timebins
%
%       smoothSignal, a value of 1 indicates that the rows of ROImtx should
%       be smoothened (moving average) before data processing, 0 indicates
%       analysis should proceed without smoothing
%

%OUTPUT: merged_ROI, new ROI signal matrix where correlated ROIs have
%        been merged and their signals have been averaged
%        roiClustIDs, information about ROIs merged to one cluster

ROIs_ts      = varargin{1,1};
smoothSignal = varargin{1,2};
sData        = varargin{1,3};
pc_rois      = varargin{1,4};

% Load deconvolved data
ROI_deconv = sData.imdata.roiSignals(2).ciaDeconvolved(pc_rois,:);


merged_ROI        = [];
merged_ROI_deconv = [];
roiClustIDs       = [];

try 
    roiClustIDs = varargin{1,5};
catch
    roiClustIDs = [];
end

% reduce noise effects
if smoothSignal == 1
    ROIs_ts_smooth = smoothdata(ROIs_ts,2,'gaussian',10);  
end

% data standardization:
ROIs_ts_zscore = zscore(ROIs_ts_smooth, [],2);

nRois       = size(ROIs_ts_zscore,1);

Y          = pdist(ROIs_ts_zscore, 'correlation');
roiVals    = squareform(Y);
temp       = tril(roiVals, -1);
uniqueCorr = temp(temp ~= 0);
    
% % cluster
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
               merged_ROI        = ROIs_ts;
               merged_ROI_zscore = ROIs_ts_zscore;
               merged_ROI_deconv = ROI_deconv;
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
                    clustSignals_zscore(t,:) = mean(ROIs_ts_zscore(rois, :), 'omitnan');
                    clustSignals(t,:)        = mean(ROIs_ts(rois, :), 'omitnan');
                    clustSignals_deconv(t,:) = mean(ROI_deconv(rois, :), 'omitnan');
                    t = t + 1;        
                end
                A = [];
                
            end
    
    %         if ~exist('clustSignals')
    %            disp('No significantly correlated clusters!')
    %            merged_ROI        = ROIs_ts;
    %            merged_ROI_zscore = ROIs_ts_zscore;
    %            merged_ROI_deconv = ROI_deconv;
    %            roiClustIDs = [];
    %         else   
                ROIs_ts(clustRois,:)        = [];
                ROIs_ts_zscore(clustRois,:) = [];
                ROI_deconv(clustRois,:)     = [];
        
                merged_ROI        = [ROIs_ts;        clustSignals];
                merged_ROI_zscore = [ROIs_ts_zscore; clustSignals_zscore];
                merged_ROI_deconv = [ROI_deconv;     clustSignals_deconv];
        
                n_unchanged_rois  = size(ROIs_ts,1);
        end

        else
               disp('No significantly correlated clusters!')
               merged_ROI        = ROIs_ts;
               merged_ROI_zscore = ROIs_ts_zscore;
               merged_ROI_deconv = ROI_deconv;
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
        clustSignals_zscore = mean(ROIs_ts_zscore(rois,:), 'omitnan');
        clustSignals        = mean(ROIs_ts(rois,:), 'omitnan');
        clustSignals_deconv = mean(ROI_deconv(rois,:), 'omitnan');
        
        % Remove clustered ROIs
        ROIs_ts(rois,:)        = [];
        ROIs_ts_zscore(rois,:) = [];
        ROI_deconv(rois,:)     = [];
    
        merged_ROI        = [ROIs_ts;        clustSignals];
        merged_ROI_zscore = [ROIs_ts_zscore; clustSignals_zscore];
        merged_ROI_deconv = [ROI_deconv;     clustSignals_deconv];
    end
end


    
    
    
        
