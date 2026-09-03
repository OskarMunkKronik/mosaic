%% Initialization

% Options.Paths.save2mat = 'N:\SCIENCE-PLEN-ECP_AnalytChem\People\Nadine\MOSAIC\260810'
[clusters_all_samples, feature_groups_all] = load_MF_results_parts( fullfile(Options.Paths.save2mat),'MF_results_part_');
feature_table = full([feature_groups_all.basePeak; clusters_all_samples'; feature_groups_all.Rt;feature_groups_all.sample;feature_groups_all.volume_summed]');
Options.processing_time.curve_resolution.process.all = tic;
X = accumarray([feature_table(:, 2),feature_table(:, 5)],feature_table(:,6));
% sum(sum(X>0,2)>=Options.do_present_in_samples)/size(X,1)
%% change
Options.do_present_in_samples = 3;
Options.blank_ratio =3;
%%
ind_keep_component_groups  = ismember(feature_table(:,2),find(sum(X>0,2) >= Options.do_present_in_samples));
feature_table = feature_table(ind_keep_component_groups,:);
clusters_all_samples = clusters_all_samples(ind_keep_component_groups);
feature_groups_all = feature_groups_all(ind_keep_component_groups);

clear X
% mz_vec = full(unique([feature_groups_all.basePeak]));
% 
% d = full([feature_groups_all.basePeak; clusters_all_samples'; ...
%     feature_groups_all.Rt; feature_groups_all.sample; ...
%     feature_groups_all.volume_summed; feature_groups_all.volume_basePeak]');
%% New initialize
Options.processing_time.curve_resolution.process.organize_data = tic;
mass_spectra_cell = {feature_groups_all.massSpec};

% Remove zeros and map mz
valid_idx = cellfun(@(s) s > 0, mass_spectra_cell(:), 'UniformOutput', false);

mz_in_mass_spectra = cellfun(@(mz,idx) mz(idx), ...
    repmat({mzroi_aug}, size(valid_idx)), valid_idx, 'UniformOutput', false);
% mz_in_mass_spectra = cellfun(@full, mz_in_mass_spectra(:), ...
%     'UniformOutput', false);

mass_spectra_cell = cellfun(@(spec,idx) spec(idx), ...
    mass_spectra_cell(:), valid_idx, 'UniformOutput', false);

mass_spectra_cell = cellfun(@full, mass_spectra_cell(:), ...
    'UniformOutput', false);
Options.Paths.save2mat = Options.Paths.save2mat;
mz_vec = full(unique(vertcat(mz_in_mass_spectra{:})));

% %% Main loop
% clc
% 
% [mz_ind, peak_borders, Rt_cluster] = deal(cell(numel(mz_vec),1));
% grps = [];
% active = true(length(clusters_all_samples),1);
% fprintf('Checking components for overlap...\n');
% q = progressParfor(numel(mz_vec));
% 
% for i = 1:numel(mz_vec)
%     % q = progressParfor(sum(active));
%     if ~active(i)
%         continue;
%     end 
%     send(q, i);
% 
%     % === Find m/z indices and borders ===
%     [mz_ind{i}, peak_borders{i},~, grps, Rt_cluster{i},active] = ...
%         get_mz_ind_for_curve_resolution([], mz_vec(i), clusters_all_samples, feature_groups_all, Options,mz_in_mass_spectra,active);
% end
%%
groups = unique(feature_table(:,2));

% mz_ind = cell(numel(groups),1);
[mz_ind, peak_borders, Rt_cluster] = deal(cell(numel(groups),1));
fprintf('Checking components for overlap...\n');

num_clusters = numel(mz_range);
q = progressParfor(num_clusters);
[x_range,y_range,mz_range] = deal(cell(num_clusters, 1));
parfor n_cluster = 1:num_clusters
        send(q, n_cluster);
    idx = clusters_all_samples == groups(n_cluster);

    % Collect unique active MS features
    % massSpecVals{n_cluster} = arrayfun(@(x) find(x.massSpec > 0), feature_groups_all(idx), 'UniformOutput', false);

    %1D
    peak_borders{n_cluster}(1,1) =  floor(median(cell2mat(arrayfun(@(x) find(x.elution1 > 0,1,'first'), feature_groups_all(idx), 'UniformOutput', false))));
    peak_borders{n_cluster}(1,2) =  ceil(median(cell2mat(arrayfun(@(x) find(x.elution1 > 0,1,'last'), feature_groups_all(idx), 'UniformOutput', false))));

    %2D
    peak_borders{n_cluster}(2,1) =  floor(median(cell2mat(arrayfun(@(x) find(x.elution2 > 0,1,'first'), feature_groups_all(idx), 'UniformOutput', false))));
    peak_borders{n_cluster}(2,2) =  ceil(median(cell2mat(arrayfun(@(x) find(x.elution2 > 0,1,'last'), feature_groups_all(idx), 'UniformOutput', false))));

    % Rt{n_cluster} = round(median(Rt_mat(:,c == n_cluster),2));
    % mz_ind{n_cluster} = unique(vertcat(massSpecVals{n_cluster}{:}));
    % Calculate unique mz indices/values for this group
    mz_range{n_cluster} = find(any(unique(vertcat(mz_in_mass_spectra{idx}))' == mzroi_aug,2));
    % mz_range{n_cluster} = unique(vertcat(mz_in_mass_spectra{idx}));
    x_range{n_cluster} = peak_borders{n_cluster}(1,1):peak_borders{n_cluster}(1,2);
    y_range{n_cluster} = peak_borders{n_cluster}(2,1):peak_borders{n_cluster}(2,2);
end

% %% Flatten mz ind
% clc
% mz_ind = flatten_cell_array(mz_ind);
% peak_borders =  flatten_cell_array(peak_borders);
% Rt_cluster =  flatten_cell_array(Rt_cluster);
% 
% %% Get ranges
% Filter = Options.Filter;
% clc
% fprintf('Processing each sample sequentially...\n');

% 
% % === Preallocate Z cell structure ===
% [x_range,y_range,mz_range] = deal(cell(num_clusters, 1));
% 
% q = progressParfor(num_clusters);
% parfor n_cluster = 1:num_clusters
%     send(q, n_cluster);
%     x_range{n_cluster} = peak_borders{n_cluster}(1,1):peak_borders{n_cluster}(1,2);
%     y_range{n_cluster} = peak_borders{n_cluster}(2,1):peak_borders{n_cluster}(2,2);
%     mz_range{n_cluster} = mz_ind{n_cluster};
% end
save(fullfile(Options.Paths.save2mat,'curve_resolution_final_meta_data_Ver2.mat'),"matFile_list","Options","mzroi_aug","x_range","y_range","mz_range","mz_ind","peak_borders","-v7.3","-nocompression");

%% get_clusters_for_deconvolution
get_clusters_for_deconvolution(matFile_list, x_range, y_range, mz_range, Options.Filter, num_clusters, n_samples, Options); % update so it takes the peak borders from MF
Options.processing_time.curve_resolution.process.organize_data = toc(Options.processing_time.curve_resolution.process.organize_data);

%% === NNMF / curve resolution ===
Options.processing_time.curve_resolution.process.nmf = tic;
q = progressParfor(num_clusters);
clc

fprintf('Curve-resolution in processes...\n');
warning('off','all')
run_cluster_chunks(matFile_list,Options, x_range, y_range, mz_range, ...
    Options.Filter, num_clusters, n_samples, 1)
warning('on','all')
Options.processing_time.curve_resolution.process.nmf = toc;

%% Build results structure
clc
Options.processing_time.curve_resolution.process.build_results_structure = tic;
[feature_groups_all, clusters_all_samples] = build_results_curve_resolution(...
    [],[] , x_range, y_range, mz_range,[], peak_borders, mz_ind, mzroi_aug, ...
    matFile_list, n_samples,Options);
feature_table = full([feature_groups_all.basePeak; clusters_all_samples; feature_groups_all.Rt;feature_groups_all.sample;feature_groups_all.volume_summed]');

Options.processing_time.curve_resolution.process.build_results_structure = toc(Options.processing_time.curve_resolution.process.build_results_structure);


%% save 
Options.processing_time.curve_resolution.save_before_duplicate = tic;
save(fullfile(Options.Paths.save2mat,'feature_table_before_duplicate_remove.mat'), "feature_table","clusters_all_samples");
save_MF_results_parts(clusters_all_samples,feature_groups_all, fullfile(Options.Paths.save2mat,'NMF_before'),'NMF_results_before_duplicate_remove_part', 10)
% save(fullfile(Options.Paths.save2mat,'feature_groups_all_before_duplicate_remove.mat'), "feature_groups_all","-v7.3","-nocompression");
Options.processing_time.curve_resolution.save_before_duplicate = toc(Options.processing_time.curve_resolution.save_before_duplicate);
%% Remove duplicates
clc
Options.processing_time.curve_resolution.process.remove_duplicates = tic;
fprintf('Removing duplicates in processes...\n');
[feature_table, feature_groups_all, clusters_all_samples,clusters_all_samples_before_duplicate_removal] = cleanup_clusters_feature_table(feature_groups_all,clusters_all_samples,feature_table, Options);

% [feature_table, feature_groups_all, clusters_all_samples,clusters_all_samples_before_duplicate_removal] = cleanup_clusters_feature_table(feature_groups_all,clusters_all_samples,feature_table, Options)
% % % % % % % % clusters_all_samples_before_duplicate_removal = clusters_all_samples;
Options.processing_time.curve_resolution.process.remove_duplicates = toc(Options.processing_time.curve_resolution.process.remove_duplicates);

%% Save results
Options.processing_time.curve_resolution.save_final = tic;
fprintf(1,'Save results\n');
save_MF_results_parts(clusters_all_samples,feature_groups_all, fullfile(Options.Paths.save2mat,'NMF_final'),'NMF_final_part',10)

% timestamp = datestr(now,'yymmdd_HHMMSS');
 % save(fullfile(Options.dir.results, ...
 %    'curve_resolution_final.mat'),...
 %     "clusters_all_samples","clusters_all_samples_before_duplicate_removal", "feature_groups_all","feature_table","-v7.3","-nocompression");

% Clear temporary variables and finalize the process
fprintf('Process completed successfully.\n');
Options.processing_time.curve_resolution.save_final = toc(Options.processing_time.curve_resolution.save_final);
Options.processing_time.curve_resolution.process.all = toc(Options.processing_time.curve_resolution.process.all);