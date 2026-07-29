%% Initialization
addpath H:\PhD\Articles\8_LCxLC_Vandalf_BigDataSet\Part_2_Data\Scripts
S_setup_load
 Options.dir.results = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob_no_MF'
close all
clc
tic
% mz_vec = full(unique([feature_groups_all.basePeak]));
mz_vec = mz_IS;
mz_vec(1) = mz_vec(1) + Options.Adducts(2) - Options.mzHydrogen;
% ;

d = full([feature_groups_all.basePeak; clusters_all_samples'; ...
    feature_groups_all.Rt; feature_groups_all.sample; ...
    feature_groups_all.volume_summed; feature_groups_all.volume_basePeak]');

%% Main loop
[mz_ind, peak_borders, Rt_cluster] = deal(cell(numel(mz_vec),1));
parfor i = 1:numel(mz_vec)

    % === Find m/z indices and borders === %
    [mz_ind{i}, peak_borders{i},~, ~, Rt_cluster{i}] = ...
        get_mz_ind_for_curve_resolution(d, mz_vec(i), clusters_all_samples, feature_groups_all, Options);

end
clear d
%% Flatten mz ind
mz_ind = flatten_cell_array(mz_ind);
peak_borders =  flatten_cell_array(peak_borders);
Rt_cluster =  flatten_cell_array(Rt_cluster);

%%
Filter = Options.Filter;
clc
fprintf('Processing each sample sequentially...\n');

% === Preallocate Z cell structure ===
num_clusters = numel(Rt_cluster);
[x_range,y_range,mz_range,mz_range_MF] = deal(cell(num_clusters, 1));

q = progressParfor(num_clusters);
parfor n_cluster = 1:num_clusters
    send(q, n_cluster);
    x_range{n_cluster} = peak_borders{n_cluster}(1,1):peak_borders{n_cluster}(1,2);
    y_range{n_cluster} = peak_borders{n_cluster}(2,1):peak_borders{n_cluster}(2,2);
    mz_range{n_cluster} = [1:length(mzroi_aug)];%mz_ind{n_cluster};
    mz_range_MF{n_cluster} = mz_ind{n_cluster};
end

%% get_clusters_for_deconvolution
 Z = get_clusters_for_deconvolution(matFile_list, x_range, y_range, mz_range, Options.Filter, num_clusters, n_samples, Options);


%% === NNMF / curve resolution ===
[W, H, cg, exp_var] = deal([]);
n_mz_vec = zeros(size(W));
q = progressParfor(num_clusters);
clc
fprintf('Curve-resolution in processes...\n');

for n_cluster = 1:num_clusters
    send(q, n_cluster);

    % Gather all Z blocks for this cluster
    Z_tmp = get_cluster_data(Options, n_cluster, n_samples, num_clusters);
   
    if ~isempty(Z_tmp)
    Z_tmp = build_4D_cluster_cube(Z_tmp, n_cluster, n_samples);
    
    %get indices
    [mz_ind{n_cluster},mz_range{n_cluster}] = deal(find(squeeze(sum(Z_tmp,[1:2,4]) > 0)));
    Z_tmp(:,:,sum(Z_tmp,[1:2,4]) == 0,:)  =[];
   
    n_mz_vec(n_cluster) = size(Z_tmp,3);
    %Curve-resolution
     [W, H, cg, exp_var,~] = perform_nnmf_unfold(Z_tmp, Options);
         % --- Refold all components once per cluster ---
            %Curve-resolution
         % --- Refold all components once per cluster ---
      % Z_refolded = refold_Z_blocks(W, H, ...
      %   x_range{n_cluster}, y_range{n_cluster}, mz_range{n_cluster}, ...
      %   n_samples, n_cluster);
    
    save_curve_resolution_results(W,H,cg, x_range{n_cluster}, y_range{n_cluster}, mz_range{n_cluster}, ...
                              n_samples, n_cluster, num_clusters, ...
                              Options,'curve_resolution_results')


    end 
 end
% Clear temporary variables and finalize the process
fprintf('Curve-resolution completed successfully.\n');



%% Build results structure
[feature_groups_all, clusters_all_samples] = build_results_curve_resolution(...
        W,H , x_range, y_range, mz_range,cg, peak_borders, mz_ind, mzroi_aug, ...
        matFile_list, n_samples,Options);
toc
%%
fprintf(1,'Save results\n');
save(fullfile(Options.dir.results,'curve_resolution_NO_MF_final.mat'), "clusters_all_samples", "feature_groups_all","-v7.3","-nocompression");

% Clear temporary variables and finalize the process
fprintf('Process completed successfully.\n');


