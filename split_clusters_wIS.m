function [newCluster,allPeaks_new] = split_clusters_wIS(mzTarget, k, allPeaks, clusters, mzroi_aug, Options,refSpec)
newCluster = clusters{k};
allPeaks_new = allPeaks{k};

% NOTE: find the real compounds m/Z and the those that are not this must be
% IS
mzTarget = [mzTarget,mzTarget(2)+Options.Adducts-Options.mzHydrogen];
% Find feature closest to target m/z
[~, fG] = min(abs(mzroi_aug(allPeaks{k}(:,4)) - mzTarget(1)));

% Get indices of all features in the same cluster
ind = clusters{k}(fG) == clusters{k};
sub_mz_cluster = allPeaks{k}(ind, 4);

% ppm window match
ind_mat = abs((mzroi_aug(sub_mz_cluster) - mzTarget) ./ mzTarget * 1e6) < Options.ppm_dev ;

if ~any(ind_mat(:,2:end))
    return
end
% Candidate group (any mz within ppm window)
% Check for insource fragment:
ind_in_source =sum(abs((full(mzroi_aug(sub_mz_cluster)-mzroi_aug(sub_mz_cluster)')- (mzTarget(1)-mzTarget(2)))./mzTarget(2)*10^6)<Options.ppm_dev,2)>0;
diff_mat = full(mzroi_aug(sub_mz_cluster)-mzroi_aug(sub_mz_cluster)')+mzTarget(2);
for nn =3:length(mzTarget)
    ind_in_source = ind_in_source+sum(abs(diff_mat-mzTarget(nn))/mzTarget(nn)*10^6<Options.ppm_dev,1)'>0;
end
if ~isempty(refSpec)
    ind_in_source = ind_in_source + sum(abs((mzroi_aug(sub_mz_cluster) - refSpec') ./ refSpec' * 1e6) < Options.ppm_dev,2)>0;
    refSpec_non_labelled = [refSpec- (mzTarget(1)-mzTarget(2));mzTarget(3:end)'];
    ind_in_source_non_labelled = sum(abs((mzroi_aug(sub_mz_cluster) - refSpec_non_labelled') ./ refSpec_non_labelled' * 1e6) < Options.ppm_dev,2)>0;

end
% Isotope set selection
% ind_IS = unique(sub_mz_cluster(ind_mat(:,1))) <= sub_mz_cluster; %Not
% % used?
% % ind_IS = false(size(ind_mat));
% ind_IS(ind_in_source,1) = true;
% % ind_IS(~ind_in_source_non_labelled,1) = true;
% ind_IS(ind_in_source,2:end) =false;
% % ind_IS(~ind_in_source_non_labelled,2:end) =false;
ind_IS = false(size(ind_mat));  % initialize all to false
ind_IS(ind_in_source, 1) = true;
ind_IS(ind_in_source, 2:end) = false;


% Map back to global indices
ind_IS_final = ismember(allPeaks{k}(:,4), sub_mz_cluster(ind_IS));
ind_adduct =  ismember(allPeaks{k}(:,4), sub_mz_cluster(~sum([ind_mat(:,2:end),ind_in_source_non_labelled],2)));

% Reassign cluster IDs for matching subset
newClusterID = max(clusters{k}) + 1;
newCluster(ind & ind_IS_final & ind_adduct ) = newClusterID;

% Flag IS
allPeaks_new(ind & ind_IS_final,7) = 1;

end

