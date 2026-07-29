clf
clc
load_dir =fullfile(Options.Paths.save2mat,'component_cluster');

u = unique(feature_table(:,2));
cluster_vec = u([1000])

%% plot peaks
cluster_vec = 3
clf
[Z_refolded,mass_spec,mz] = deal(cell(max(cluster_vec),1));
figure(1)
for n = 1:length(cluster_vec)
    n_cluster = cluster_vec(n);
    [~, ~, Z_refolded(n_cluster), mass_spec(n_cluster), mz(n_cluster), p_tmp] = plotClusters_from_disk( ...
        n_cluster, feature_table, clusters_all_samples, ...
        mzroi_aug, load_dir, Options, n_samples, ...
        max(rand([110,3]),[0,0,0]), 1:n_samples, false);
end 

%% Mass spec
figure(2)
stem(mz{n_cluster},mass_spec{n_cluster},'Marker','none','LineWidth',1)
xlabel('m/z')
ylabel('Intensity')