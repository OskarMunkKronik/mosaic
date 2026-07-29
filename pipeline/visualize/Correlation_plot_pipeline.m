
R = corr(T', 'Rows', 'pairwise');

% Replace NaNs (from constant columns) with zero correlation
R(isnan(R)) = 0;

% --- Convert correlation to distance ---
% Use correlation distance: d = 1 - |r|
D = 1 - abs(R);

% Force diagonal = 0
D(1:size(D,1)+1:end) = 0;

% --- Convert to condensed pdist vector (REQUIRED by linkage & optimalleaforder) ---
D_vec = squareform(D);   % condensed form

% --- Perform hierarchical clustering ---
Z = linkage(D_vec, 'average');

% --- Optimize leaf order ---
order = optimalleaforder(Z, D_vec);

% --- Reorder correlation matrix ---
R_sorted = R(order, order);

% subplot(1,3,1)
clf
imagesc(R_sorted);
set(gca,'YDir','normal')
colormap(turbo);
colorbar;

%
