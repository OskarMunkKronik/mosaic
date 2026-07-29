mzTarget  =  314.1391

mz_IS(14);
k = 2;

clf
col = lines(78);
EIC = cell(78,1);
%%

for k = 1:78
    n = get_cluster_members(feature_groups_all, clusters_all_samples, mzTarget, k,Options);

    if ~isempty(n)
        subplot(2,1,1)
        [h,EIC{k}] = plotFeatureGroup3D(feature_groups_all, Options, n,'color',col(k,:));
        vol(k) = feature_groups_all(n).volume_summed;

        hold on
        subplot(2,1,2)
        h = plotMassSpec(feature_groups_all, mzroi_aug, n,1,false, 'color',col(k,:));
        hold on
    end
end
