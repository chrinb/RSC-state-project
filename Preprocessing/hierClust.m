function [merged_ROI,roiClustIDs, interClust, intraClust,  tree, uniqueCorr ] = hierClust(file, ROImtx, smoothSignal, roiArray)

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

%smoothen the ROI signals
% minSpeed   = settings.minspeed;
% maxSpeed   = settings.maxspeed; %*max(sData.behavior.runSpeedDs(:));
file = strcat(file, '/clustering/');
if ~exist(file)
    mkdir(file)
end

linkg      = 'average'; %settings.linkage.Selection;

% reduce noise effects
ROIsignals = zeros(size(ROImtx));
if smoothSignal == 1
    ROIsignals = smoothdata(ROImtx,2,'gaussian',10);
else ROIsignals = ROImtx;
end

% data standardization:
ROIsignals = zscore(ROIsignals, [],2);

merged_ROI = []; roiClustIDs = []; interClust = []; intraClust = []; 
nRois = size(ROIsignals,1);

% waitbar(0.1,wait,'Preprocessing ist finished. Calculation of correlation matrix...')
%roiVals = corrcoef(ROIsignals','Rows','pairwise'); 
% waitbar(0.9,wait,'Calculation of correlation matrix is finished. Clustering...')

Y = pdist(ROIsignals, 'correlation');
roiVals = squareform(Y);
temp = tril(roiVals, -1);
uniqueCorr = temp(temp ~= 0);

    
% % cluster
tree = linkage(double(Y), linkg);%, 'correlation'); 
for i = 1:size(ROIsignals,1)
     for j = i+1:size(ROIsignals,1)
         if i == j      
         else
                 centerAll(1,:) = roiArray(1, i).center;
                 centerAll(2,:) = roiArray(1, j).center;
                 distValTemp(i,j)   = pdist(centerAll, 'euclidean');
         end
     end
end


% figure()
% H= dendrogram(tree, size(tree,1));
% %dendrogram(Z, size(Z,1), 'Orient', 'Left', 'Labels', species);
% set(H,'LineWidth',1.5, 'Color', [0 0 0])
% yline(0.3, 'LineWidth', 2, 'LineStyle','--')
% box on
% set(gca, 'FontSize', 20, 'FontName', 'Gotham')
% set(gca, 'XTickLabel', [])
% ylabel('Distance')
% saveas(gca,strcat(file, 'dendrogram.fig'),'fig')
% close


T = cluster(tree,'Cutoff', 0.3, 'Criterion', 'distance');
nclust = max(T);

% sort the clusters by highest mean correlation
 nclustReal = []; j = 1; meanClust = [];
for iclust=1:nclust   
    % calculate intra cluster similarity
    if length(find(T== iclust)) > 1
        vla             = nanmean(roiVals(T==iclust,T==iclust));
        nreal(j)        = nanmean(nanmean(roiVals(T==iclust,T==iclust)));
        meanClust(j,:)  = nanmean(ROIsignals(T==iclust, :));
        indices         = find(T==iclust);

        for i = 1:length(indices)
            vals = [meanClust(j,:); ROIsignals(indices(i),:)];
            distVal(i) = pdist(vals, 'correlation'); 
        end
        
        intraClust(j)      = nanmean(distVal);
        intraClustStd(j)   = nanstd(distVal);
        distVal = [];
        nclustReal = [nclustReal iclust];
        j = j+1;
    end
end

% calculate inter cluster similarity
if isempty(meanClust) || size(meanClust,1) == 1 ; return; end
meanDistClust = squareform(pdist(meanClust, 'correlation'));
for i = 1:length(nclustReal(:))
    interClust(i) = nanmean(meanDistClust(i,:));
    interClustStd(i) = nanstd(meanDistClust(i,:));
end

[qsort_clust,indclust] = sort(nreal,'descend');

nclust = length(nclustReal);
indcl = [];
for iclust=1:nclust
    cluster_pre = transpose(find(T==indclust(iclust)));

    for i_entry = 1:length(cluster_pre)
        unsort_entry(i_entry) = mean(roiVals(i_entry,:));%roiCorr_raw(i_entry,:));
    end
    [qsort_entry,indclust_entry] = sort(unsort_entry,'descend');
    cluster_post = cluster_pre(indclust_entry);
        
    indcl = [indcl cluster_post];
    indclust_entry =[];
    cluster_post   =[];
    cluster_pre    =[];
    unsort_entry   =[];
end

finalClust = indclust(qsort_clust < 0.3);%corrThresh);
finalClust = nclustReal(finalClust);

mean(qsort_clust ~= 1)

% merge axons within highest correlated clusters
T(:,2) = 1:nRois;
T = sortrows(T,1);
clustRois = [];
if isempty(finalClust)
    disp('No significantly correlated clusters!')
    merged_ROI = [];
else
    t = 1;
    for i = 1:length(finalClust)
 
        A = find(T(:,1) == finalClust(i));
        if length(A) > 1
            clustRois                      = [clustRois; T(A,2)];
            roiClustIDs(i).rois            = T(A,2);            
            clustSignals(t,:) = nanmean(ROIsignals(T(A,2),:));
            t = t + 1;        
        end
        A = [];
        
    end
    if ~exist('clustSignals')
       disp('No significantly correlated clusters!')
       merged_ROI = ROIsignals;
       roiClustIDs = [];
    else   
        ROIsignals(clustRois,:) = [];
        merged_ROI = [ROIsignals; clustSignals];
    end
end
    
