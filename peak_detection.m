addpath C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection
mzTarget = 268.191;
[~, a] = min(abs(mzroi_aug - mzTarget));
Options.Rtdev(1:2) = [5,15];

%%
% Initialize results as empty matrix
allPeaks = [];
k = 3
[peak_borders,sum_elution_profile,Z_feature]  = deal(cell(length(mzroi_aug),1));
close all
tic
sz = size(double(Z(:,:,1)));
n_feature = 1;
a = find(max(MSroi_aug{k},[],2)>Options.int_thresh);
for nn = a-10:a+100;%1:length(a)%a-1000:a+1000
    n = a(nn);
    % Extract the n-th slice (3rd mode)
    data = double(Z(:,:,n));

    if any(data, "all")
        [pks, locs_y, locs_x] = peaks2((smoothdata(data,2,'movmean',5)'), 'MinPeakHeight', Options.int_thresh,'MinPeakDistance',3);

        if ~isempty(pks)
            for n_peak = 1:length(locs_x)
                [peak_borders{n_feature},sum_elution_profile_tmp,Z_feature{n_feature}] = get_rt_square(locs_x(n_peak), locs_y(n_peak), data, Options);
                for dim =1:2
                    sum_elution_profile{n_feature}{dim} = zeros(sz(dim),1);
                    sum_elution_profile{n_feature}{dim}(peak_borders{n_feature}{dim}(1):peak_borders{n_feature}{dim}(2)) = sum_elution_profile_tmp{dim};

                end
                %             [peak_borders,sum_elution_profile,Z_feature]
                %                             Z_feature
                n_feature = n_feature+1;
            end
            % Concatenate results into a matrix
            temp = [pks(:), locs_y(:), locs_x(:), repmat(n, numel(pks), 1)];
            allPeaks = [allPeaks; temp];
        end
    end
end
n_feature = n_feature-1;
peak_borders = peak_borders(1:n_feature);
sum_elution_profile = sum_elution_profile(1:n_feature);
Z_feature = Z_feature(1:n_feature) ;


%%
nFeat = numel(sum_elution_profile);

% Preallocate similarity matrices
cosMat = cell(1,2);
for dim = 1:2
    cosMat{dim} = zeros(nFeat); % square similarity matrix
end
cosMat_combined = zeros(nFeat);

% Compute cosine similarities
for nn = 1:nFeat
    for n = 1:nFeat
        for dim = 1:2
            x = sum_elution_profile{nn}{dim}(:)./sum(sum_elution_profile{nn}{dim},'all');
            y = sum_elution_profile{n}{dim}(:)./sum(sum_elution_profile{nn}{dim},'all');
            cosMat{dim}(n,nn) = dot(x,y) / (norm(x)*norm(y) + eps);
        end
        cosMat_combined(n,nn) = cosMat{1}(n,nn) * cosMat{2}(n,nn);
    end
end

% --- Convert similarity to distance ---
distMat = 1 - cosMat_combined;
distMat = (distMat + distMat')/2;                % enforce symmetry
distMat(1:nFeat+1:end) = 0;                      % zero diagonal
distVec = squareform(distMat);                   % condensed form

% --- Hierarchical clustering ---
Z_clust = linkage(distVec, 'average');
threshold = 0.3;
clusters = cluster(Z_clust, 'cutoff', threshold, 'criterion', 'distance');

% --- Display results ---
numClusters = max(clusters);
disp(['Number of clusters: ', num2str(numClusters)]);

figure;
% cosMat_combined(cosMat_combined==0) = nan;
imagesc(cosMat_combined);
% colormap hot
colorbar ;
axis square;
title('Combined Cosine Similarity');
xlabel('Features');
ylabel('Features');

figure;
gscatter(mzroi_aug(allPeaks(:,4)), clusters, clusters); 
xlabel('m/z');
ylabel('Cluster ID');
title('Cluster Assignment by m/z');

%%
close all
mzTarget =268.191
[~,a] = min(abs(mzroi_aug(allPeaks(:,4))-mzTarget));
feature_vec =     find(clusters(a) == clusters);

mass_spectrum = zeros(length(feature_vec),1);
subplot(1,2,1)
for nn = 1:length(feature_vec)
    n = feature_vec(nn);
        [X_coord,Y_coord] = meshgrid(peak_borders{n}{1}(1):peak_borders{n}{1}(2), ...
            peak_borders{n}{2}(1):peak_borders{n}{2}(2));
%     Z_recon  = sum_elution_profile{n}{1}*sum_elution_profile{n}{2}';
mass_spectrum(nn) = sum( Z_feature{n},'all');
contour(X_coord,Y_coord,        Z_feature{n}')
hold on
end 

subplot(1,2,2)
stem(mzroi_aug(allPeaks(feature_vec,4)),mass_spectrum)

%%


