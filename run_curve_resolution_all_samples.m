%% Initialization
addpath H:\PhD\Articles\8_LCxLC_Vandalf_BigDataSet\Part_2_Data\Scripts
S_setup_load
close all
clc
tic;
mz_vec = full(unique([feature_groups_all.basePeak]));
% mz_vec(1) = mz_vec(1) + Options.Adducts(2) - Options.mzHydrogen;

d = full([feature_groups_all.basePeak; clusters_all_samples'; ...
    feature_groups_all.Rt; feature_groups_all.sample; ...
    feature_groups_all.volume_summed; feature_groups_all.volume_basePeak]');

Options.dir.results = fullfile('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob','nnmf_251119_ver3');

%% Main loop
clc
[mz_ind, peak_borders, Rt_cluster] = deal(cell(numel(mz_vec),1));
grps = [];
active = true(length(clusters_all_samples),1);
fprintf('Checking components for overlap...\n');
q = progressParfor(numel(mz_vec));
% mass_spectrum_tmp = arrayfun(@(x) find(x.massSpec > 0), feature_groups_all, 'UniformOutput', false);
parfor i = 1:numel(mz_vec)

    send(q, i);
    % if active(i)
    % === Find m/z indices and borders ===
    [mz_ind{i}, peak_borders{i},~, grps, Rt_cluster{i}] = ...
        get_mz_ind_for_curve_resolution(d, mz_vec(i), clusters_all_samples, feature_groups_all, Options,mass_spectrum_tmp,mzroi_aug);
    % disable = ismember(clusters_all_samples,grps);
    % active(disable) = false;
end
% end
%% Flatten mz ind
clc
mz_ind = flatten_cell_array(mz_ind);
peak_borders =  flatten_cell_array(peak_borders);
Rt_cluster =  flatten_cell_array(Rt_cluster);

%% Get ranges
Filter = Options.Filter;
clc
fprintf('Processing each sample sequentially...\n');
num_clusters = numel(mz_ind);
% === Preallocate Z cell structure ===
[x_range,y_range,mz_range] = deal(cell(num_clusters, 1));

q = progressParfor(num_clusters);
parfor n_cluster = 1:num_clusters
    send(q, n_cluster);
    x_range{n_cluster} = peak_borders{n_cluster}(1,1):peak_borders{n_cluster}(1,2);
    y_range{n_cluster} = peak_borders{n_cluster}(2,1):peak_borders{n_cluster}(2,2);
    mz_range{n_cluster} = mz_ind{n_cluster};
end
save(fullfile(Options.dir.results,'curve_resolution_final_meta_data_Ver2.mat'),"matFile_list","Options","mzroi_aug","x_range","y_range","mz_range","mz_ind","peak_borders","-v7.3","-nocompression");

%% get_clusters_for_deconvolution
clc
Z = get_clusters_for_deconvolution(matFile_list, x_range, y_range, mz_range, Options.Filter, num_clusters, n_samples, Options);

%% === NNMF / curve resolution ===
[W, H, cg, exp_var] = deal([]);
q = progressParfor(num_clusters);
clc
fprintf('Curve-resolution in processes...\n');

parfor n_cluster = 1:num_clusters
    %Progress bar
    send(q, n_cluster);

    % Gather all Z blocks for this cluster
    Z_tmp = get_cluster_data(Options, n_cluster, n_samples, num_clusters);

    if ~isempty(Z_tmp)
        Z_tmp = build_4D_cluster_cube(Z_tmp, n_cluster, n_samples);

        %Curve-resolution
        [W, H, cg, exp_var,~] = perform_nnmf_unfold(Z_tmp, Options);

        % --- Refold all components once per cluster ---
        save_curve_resolution_results(W,H,cg, x_range{n_cluster}, y_range{n_cluster}, mz_range{n_cluster}, ...
            n_samples, n_cluster, num_clusters, ...
            Options,'curve_resolution_results')

    end
end




%% Build results structure
clc
[feature_groups_all, clusters_all_samples] = build_results_curve_resolution(...
    [],[] , x_range, y_range, mz_range,[], peak_borders, mz_ind, mzroi_aug, ...
    matFile_list, n_samples,Options);
feature_table = full([feature_groups_all.basePeak; clusters_all_samples; feature_groups_all.Rt;feature_groups_all.sample;feature_groups_all.volume_summed]');
% save(fullfile(Options.dir.results,'feature_table_before_duplicate_remove.mat'), "feature_table","-v7.3","-nocompression");
%% Remove duplicates
clc
% fprintf('Extract elu_profiles...\n');
% num_clusters = max(feature_table(:,2));
% [elu1,elu2,mass_spec_idx] = deal(cell(num_clusters,1));
% load_dir = fullfile(Options.dir.results,'component_cluster');
% 
% q = progressParfor(num_clusters);
% numFiles  = numel(dir([load_dir,'\*.mat']));
% n_digits_clusters = length(num2str(numFiles));
% parfor n_cluster = 1:num_clusters
%     %Progress bar
%     send(q, n_cluster);
% 
%     [~,~,elu1{n_cluster},elu2{n_cluster},mass_spec_idx{n_cluster}] = get_SEP_massSpecIdx(n_cluster,load_dir,Options,n_digits_clusters);
% 
% end
% %%
fprintf('Removing duplicates in processes...\n');
[feature_table, feature_groups_all, clusters_all_samples,clusters_all_samples_before_duplicate_removal] = cleanup_clusters_feature_table(feature_groups_all,clusters_all_samples,feature_table,mass_spec_idx, Options);
process_time.curve_resolution = toc;

%% Save results
tic
fprintf(1,'Save results\n');
save(fullfile(Options.dir.results,'curve_resolution_final_251110_Ver2.mat'), "clusters_all_samples","clusters_all_samples_before_duplicate_removal", "feature_groups_all","feature_table","-v7.3","-nocompression");
% save(fullfile(Options.dir.results,'curve_resolution_final_meta_data_Ver2.mat'),"matFile_list","Options","mzroi_aug","x_range","y_range","mz_range","mz_ind","peak_borders","-v7.3","-nocompression");

% Clear temporary variables and finalize the process
fprintf('Process completed successfully.\n');

process_time.save_resultes = toc;