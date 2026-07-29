%% Initialization
tic;
% mz_vec = full(unique([feature_groups_all.basePeak]));
% 
% d = full([feature_groups_all.basePeak; clusters_all_samples'; ...
%     feature_groups_all.Rt; feature_groups_all.sample; ...
%     feature_groups_all.volume_summed; feature_groups_all.volume_basePeak]');
%% New initialize
mass_spectra_cell = {feature_groups_all.massSpec};

% Remove zeros and map mz
valid_idx = cellfun(@(s) s > 0, mass_spectra_cell, 'UniformOutput', false);

mz_in_mass_spectra = cellfun(@(mz,idx) mz(idx), ...
    repmat({mzroi_aug}, size(valid_idx)), valid_idx, 'UniformOutput', false);

mass_spectra_cell = cellfun(@(spec,idx) spec(idx), ...
    mass_spectra_cell, valid_idx, 'UniformOutput', false);

Options.dir.results = Options.Paths.save2mat;
mz_vec = full(unique(vertcat(mz_in_mass_spectra{:})));
%% Main loop
clc
[mz_ind, peak_borders, Rt_cluster] = deal(cell(numel(mz_vec),1));
grps = [];
active = true(length(clusters_all_samples),1);
fprintf('Checking components for overlap...\n');
q = progressParfor(numel(mz_vec));

for i = 1:numel(mz_vec)
    if ~active(i)
        continue;
    end 
    send(q, i);
    
    % === Find m/z indices and borders ===
    [mz_ind{i}, peak_borders{i},~, grps, Rt_cluster{i},active] = ...
        get_mz_ind_for_curve_resolution([], mz_vec(i), clusters_all_samples, feature_groups_all, Options,mz_in_mass_spectra,active);
end

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
get_clusters_for_deconvolution(matFile_list, x_range, y_range, mz_range, Options.Filter, num_clusters, n_samples, Options); % update so it takes the peak borders from MF

%% === NNMF / curve resolution ===
q = progressParfor(num_clusters);
clc

fprintf('Curve-resolution in processes...\n');
warning('off','all')
run_cluster_chunks(Options, x_range, y_range, mz_range, ...
    Filter, num_clusters, n_samples, 1)
warning('on','all')

%% Build results structure
clc
[feature_groups_all, clusters_all_samples] = build_results_curve_resolution(...
    [],[] , x_range, y_range, mz_range,[], peak_borders, mz_ind, mzroi_aug, ...
    matFile_list, n_samples,Options);
feature_table = full([feature_groups_all.basePeak; clusters_all_samples; feature_groups_all.Rt;feature_groups_all.sample;feature_groups_all.volume_summed]');
save(fullfile(Options.dir.results,'feature_table_before_duplicate_remove.mat'), "feature_table","-v7.3","-nocompression");
save(fullfile(Options.dir.results,'feature_groups_all_before_duplicate_remove.mat'), "feature_groups_all","-v7.3","-nocompression");

%% Remove duplicates
clc
fprintf('Removing duplicates in processes...\n');
[feature_table, feature_groups_all, clusters_all_samples,clusters_all_samples_before_duplicate_removal] = cleanup_clusters_feature_table(feature_groups_all,clusters_all_samples,feature_table, Options);
process_time.curve_resolution = toc;

%% Save results
tic
fprintf(1,'Save results\n');
% timestamp = datestr(now,'yymmdd_HHMMSS');
 save(fullfile(Options.dir.results, ...
    'curve_resolution_final.mat'),...
     "clusters_all_samples","clusters_all_samples_before_duplicate_removal", "feature_groups_all","feature_table","-v7.3","-nocompression");

% Clear temporary variables and finalize the process
fprintf('Process completed successfully.\n');

process_time.save_resultes = toc;