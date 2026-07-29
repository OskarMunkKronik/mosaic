function clusters = run_clustering(cosMat_combined,threshold)
% --- Convert similarity to distance ---
distMat = 1 - cosMat_combined;
distMat = (distMat + distMat')/2;                % enforce symmetry
distMat(1:size(distMat,1)+1:end) = 0;                      % zero diagonal
distVec = squareform(distMat);                   % condensed form

% --- Hierarchical clustering ---
Z_clust = linkage(distVec, 'average');
clusters = cluster(Z_clust, 'cutoff', threshold, 'criterion', 'distance');
end
