function [Labels]=hierachical_clastering_PG(CorrMatrix,labels_for,CustomColormap)

% Step 1: Minkowski exponent -> enhance distances
p = 9.5; %for rat10 9.5

distMatrix = (1 - (CorrMatrix)).^p; % Convert correlation to distance (1 - correlation)

% Step 2: Compute the distance matrix for columns
distMatrix=distMatrix-diag(diag(distMatrix));
distVector = squareform(distMatrix, 'tovector'); % Extract lower triangle for linkage

%Step 3: Perform hierarchical clustering to get the linkage matrix
linkageMatrix = linkage(distVector, 'ward'); % 'average' linkage is common

% Step 4: Compute the optimal threshold
optimalThreshold = findOptimalThresholdElbow(linkageMatrix);

% Visualize as a dendrogram
dendrOrder=optimalleaforder(linkageMatrix,distVector)
figure;
subplot(1,2,1)
[hDendrogram,~,order] =dendrogram(linkageMatrix,0,'ColorThreshold',optimalThreshold,'Reorder',flip(dendrOrder));
xticks(1:length(order))
xticklabels(labels_for(order))
view([-90 90])
yline(optimalThreshold,LineStyle="--",LineWidth=2)
title('Hierarchical Clustering', 'of Neuronal Spike-Count Correlations');
axis square
xtickangle(45);


subplot2=subplot(1,2,2)

plot_corr_matrix=(CorrMatrix - diag(diag(CorrMatrix)));

imagesc(plot_corr_matrix(dendrOrder,dendrOrder));
xlim([0.5 length(order)+0.5])
ylim([0.5 length(order)+0.5])
clim([-1 1])
yticks(1:length(order));
yticklabels((labels_for(dendrOrder)));
colormap(CustomColormap)
axis square
title('Pairwise Neuronal Spike-Count Correlation', '(reordered)')
colorbar(subplot2,'Position',...
    [0.915865856651938 0.357142857142857 0.0238392857142855 0.317857142857143]);
clim([-0.3 0.3])
fig=gcf;
fig.Renderer= "painters"
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 3840 2160];
fig.PaperSize = [3840 2160];


% Step 5: Identify Leaf Nodes and Extract Their Colors

% Number of leaves (original data points)
numLeaves = size(linkageMatrix, 1) + 1;

% Initialize cluster labels
clusterLabels = zeros(numLeaves, 1);

% Extract the color of each branch at the leaf level
leafColors = nan(numLeaves, 3); % Store RGB color for each leaf

for i = 1:length(hDendrogram)
    % Get branch coordinates
    XData = get(hDendrogram(i), 'XData');
    YData = get(hDendrogram(i), 'YData');
    BranchColor = get(hDendrogram(i), 'Color'); % Get branch color

    % Identify leaf nodes (where YData == 0)
    leafIndices = find(YData == 0);

    for j = 1:length(leafIndices)
        leafPos = round(XData(leafIndices(j))); % Approximate leaf index
        if leafPos >= 1 && leafPos <= numLeaves
            leafColors(leafPos, :) = BranchColor; % Store RGB color
        end
    end
end

% Find cluster labels to match the dendrogram leaf
[~, ~, clusterLabels] = unique(leafColors, 'rows', 'stable');

% Save as an output the extracted cluster labels

Labels=labels_for(dendrOrder)';
Labels(:,2)=num2cell(flip(clusterLabels));

end

%Step 4: Define a function to compute the optimal threshold (Elbow Method)
function optimalThreshold = findOptimalThresholdElbow(linkageMatrix)
    distances = linkageMatrix(:, 3); % Extract distances (third column of linkage matrix)
    deltas = diff(distances); % Compute first-order differences
    [~, elbowIndex] = max(deltas); % Find the largest difference (elbow point)
    optimalThreshold = distances(elbowIndex); % Threshold at the elbow poin

end

% 
% function optimalThreshold = findOptimalThresholdElbow(linkageMatrix)
% 
% distances = linkageMatrix(:, 3);
% % The Kneedle algorithm is a method for detecting the “elbow” or “knee” point in a monotonically increasing (or decreasing) curve,
% % where the rate of change sharply shifts. It works by comparing the original curve to a reference line (typically y = x) and finding the point of maximum deviation,
% % which corresponds to the area of highest curvature.
% % Normalize distances to [0, 1]
% x = (1:length(distances))';
% xNorm = (x - min(x)) / (max(x) - min(x));
% yNorm = (distances - min(distances)) / (max(distances) - min(distances));
% 
% % Ideal line: y = x (straight diagonal)
% idealLine = xNorm;
% 
% % Find difference between curve and diagonal
% diffFromIdeal = yNorm - idealLine;
% 
% % Find point of maximum difference (elbow) (inverted)
% [~, elbowIdx] = min(diffFromIdeal);
% plot( x ,diffFromIdeal)
% % Output the corresponding threshold
% optimalThreshold = distances(elbowIdx);
% end