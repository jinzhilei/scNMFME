function K=gap_cluster(data)
% gap_cluster Compute the optimatical cluster number in single cell matrix
% INPUT:
%   data: the single cell matrix
% OUTPUT:
%   K: the optimatical cluster number
% PCA and visualization
% [~,score] = pca(zscore(data)');
[~,score] = pca(zscore(data));
score=score(:,1:2);
eva = evalclusters(score,@kmedoids,'gap','KList',1:20,'Distance','correlation');
plot(eva)
K = eva.OptimalK;
