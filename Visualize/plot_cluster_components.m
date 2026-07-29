% tmp = load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\component_cluster\results_02388.mat")
load_dir = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\component_cluster\';
addpath C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection
load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\curve_resolution_final_251110_Ver2.mat",...
    "clusters_all_samples_before_duplicate_removal","feature_table");
tmp = load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\feature_table_before_duplicate_remove.mat")

% [ ~, ~, ~,clusters_all_samples_before_duplicate_removal] = cleanup_clusters_feature_table([],tmp.feature_table(:,2),tmp.feature_table, Options);
% save("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\curve_resolution_final_251110_Ver2.mat","clusters_all_samples_before_duplicate_removal","-append")
%%
%final cluster 
clc
for n = 1:length(clusters_09)
    if any(clusters_09(n).members == 527)
        n
    end 
end 
%%


n_cluster_vec = clusters_05(27).names%clusters.members; final_features(order([1256:1271]));
figure(4)
clf
h = 1;
clc
col = lines(length(n_cluster_vec));
[X,Y,Z_refolded,mz_vals_cluster,mass_spec] = deal(cell(length(n_cluster_vec),1)) ;
for n = 1:length(n_cluster_vec)
    subplot(2,1,h)


% Cluster numbers before duplicate removal
ind_cluster = feature_table(:,2) == n_cluster_vec(n);

%Find the cluster with most counts
[gc,grps] = groupcounts(clusters_all_samples_before_duplicate_removal(ind_cluster));


n_cluster_old = grps(1);
cd(load_dir)
numFiles  = numel(dir('*.mat'));
n_digits_clusters = length(num2str(numFiles));

% Zero-padded cluster folder name
    cluster_str = sprintf(['%0', num2str(n_digits_clusters), 'd'], n_cluster_old);
    % Output filename (zero-padded and structured)
    filename = sprintf("results_%s.mat", cluster_str);

    % === Save ===
    results = load(fullfile(load_dir, filename));
    % --- Setup indices ---
    [X{n},Y{n}] = meshgrid(Options.coords.X(1,results.x_range),Options.coords.Y(results.y_range,1));

    % [X,Y] = meshgrid(results.x_range,results.y_range);
   



    % --- Refold all components once per cluster ---
    Z_refolded{n} = refold_Z_blocks(results.W, results.H, ...
        results.x_range, results.y_range, results.mz_range, ...
        n_samples, []);
    %
    mz_vals_cluster{n} = mzroi_aug(results.mz_range);
    % clf
    % figure(4)
    for k = 1:size(Z_refolded{n},4)
        
    plot3(X{n},Y{n},squeeze(max(Z_refolded{n}(:,:,:,k),[],[3])),'Color',col(n,:))
    hold on 
    % pause
    mass_spec{n}(:,k) = squeeze(sum(Z_refolded{n}(:,:,:,k),[1:2]));
    end 
    % h = h +1;
    subplot(2,1,2)
    % squeeze(sum(Z_refolded{n},[1:2,4]))
    stem(mz_vals_cluster{n},median(mass_spec{n},2),'Color',col(n,:))
    hold on
    % h = h +1;
    % xlim([255 270])
end 
    
%%
clc
q = progressParfor(num_clusters);
clc
fprintf('Curve-resolution in processes...\n');
cd(load_dir)
numFiles  = numel(dir('*.mat'));
n_digits_clusters = length(num2str(numFiles));

parfor n_cluster = 1:max(feature_table(:,2))
    %Progress bar
    send(q, n_cluster);    
ind_cluster = feature_table(:,2) == n_cluster;

%Find the cluster with most counts
[gc,grps] = groupcounts(clusters_all_samples_before_duplicate_removal(ind_cluster));


n_cluster_old = grps(1);

% Zero-padded cluster folder name
    cluster_str = sprintf(['%0', num2str(n_digits_clusters), 'd'], n_cluster_old);
    % Output filename (zero-padded and structured)
    filename = sprintf("results_%s.mat", cluster_str);

    % === Load ===
  mz_range = load(fullfile(load_dir, filename),'mz_range');
num_mz(n_cluster) = length(mz_range.mz_range);
end 

save(fullfile(Options.dir.results,'num_fragments_in_cluster.mat'),'num_mz')
%%
[~,idx] = sort(num_mz);
int = accumarray(feature_table(:,2),feature_table(:,6),[],@median);


clf
scatter(int(idx),num_mz(idx))
xline(1000)
yline(1.5)