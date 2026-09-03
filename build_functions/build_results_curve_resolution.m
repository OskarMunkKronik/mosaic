function [Results, cluster_curve_resolution] = build_results_curve_resolution( ...
    W, H, x_range, y_range, mz_range, cg, peak_borders, mz_ind, mzroi_aug, ...
    matFile_list, n_samples,Options)

fprintf('Building Results structure...\n');

num_clusters = numel(mz_range);
load_dir_curve = fullfile(Options.Paths.save2mat, 'curve_resolution_results');
% load_dir_MF = fullfile(Options.dir.results, 'clusters_for_curve_resolution');
n_digits_clusters = ceil(log10(num_clusters ) + 1);

% === Preallocate Results ===
emptyResult = struct( ...
    'Rt', [], ...
    'elution1', [], ...
    'elution2', [], ...
    'feature_elu_2D', [], ...
    'massSpec', [], ...
    'volume_summed', [], ...
    'volume_basePeak', [], ...
    'basePeak', [], ...
    'sample_name', [], ...
    'sample', [], ...
    'sample_feature_ID', [], ...
    'sample_cluster_ID', [], ...
    'IS_flag', [], ...
    'num_components', []);
Results(num_clusters,1) = emptyResult;

cluster_curve_resolution = zeros(1, num_clusters);
counter = 1;
max_cluster = 0;
q = progressParfor(num_clusters);
% === Loop over clusters ===
for n_cluster =  1:num_clusters
    % n_cluster=vec(n_h);
    send(q, n_cluster);
    
    % Zero-padded cluster folder name
    cluster_str = sprintf(['%0', num2str(n_digits_clusters), 'd'], n_cluster);
    % Output filename (zero-padded and structured)
    filename = sprintf("results_%s.mat", cluster_str);


    if ~exist(fullfile(load_dir_curve, filename), "file")
        continue
    end

    % === Save ===
    % load(fullfile(load_dir, filename), ...
    %     "W","H","cg");
     load(fullfile(load_dir_curve, filename), "cg","Z_refolded");

    % --- Setup indices ---
    n_vec = unique(cg);
    F = max(n_vec);
    X = peak_borders{n_cluster}(1,1):peak_borders{n_cluster}(1,2);
    Y = peak_borders{n_cluster}(2,1):peak_borders{n_cluster}(2,2);



    % --- Refold all components once per cluster ---
    % Z_refolded = refold_Z_blocks(W, H, ...
    %     x_range{n_cluster}, y_range{n_cluster}, mz_range{n_cluster}, ...
    %     n_samples, n_cluster);
    %
    mz_vals_cluster = mzroi_aug(mz_range{n_cluster});
    len_mzroi =length(mzroi_aug);
    % === Loop over each subcomponent ===
    for n = n_vec(:)'  % ensures row vector
        cluster_mask = (cg == n);

        % Collapse x,y dims once — vectorized
        tmp = squeeze(sum(Z_refolded(:,:,cluster_mask,:), [1 2]));  % (mz × sample)
        mz_ind_tmp = (mz_range{n_cluster}(cluster_mask)) ;
        mass_spec = reshape(tmp, [], n_samples);

        % --- Base peak for each sample ---
        % [~, basepeak_idx] = max(mass_spec, [], 1);
        basePeak_mz = mz_vals_cluster(cluster_mask);
        % basePeak_mz = basePeak_mz(basepeak_idx);
        [~,basepeak_idx] = max(sum(mass_spec,2));
        basePeak_mz = basePeak_mz(basepeak_idx);
        component_number =n + max_cluster;
        BPC = squeeze(sum(Z_refolded(:,:,basepeak_idx,:),4));
        % save_curve_resolution_results(W,H(:,cluster_mask),n, X, Y, mz_ind_tmp, ...
        %     n_samples, component_number, num_clusters, ...
        %     Options,'component_cluster') % Save H an W as sparse matrices
           save_curve_resolution_results([],[],n,Z_refolded(:,:,cluster_mask,:), X, Y, mz_ind_tmp,n_samples, component_number, num_clusters,Options,'component_cluster') % Save H an W as sparse matrices
        
        %%
        % --- Extract base peak chromatograms (vectorized-ish) ---
        vol_tmp = squeeze(sum(Z_refolded(:,:,cluster_mask,:), 1:3));
        [r, c] = find(BPC == max(BPC(:)), 1);
        for k = 1:n_samples % build feature table directly!  this is suboptimal make this vectorized and parrallel! 
            % bp_idx = basepeak_idx(k);
            % BPC = Z_refolded(:,:,bp_idx,k);
            
            % if ~any(BPC(:)), continue; end

            % Max position
            
        
            % Fill Results
            Results(counter).Rt                          = [X(r); Y(c)];
            Results(counter).elution1                    = X;
            Results(counter).elution2                    = Y;
            % Results(counter).feature_elu_2D              = sptensor(BPC);
            Results(counter).massSpec                    = sparse(len_mzroi,1);
            Results(counter).massSpec(mz_ind_tmp)        = mass_spec(:,k);
            Results(counter).volume_summed               = vol_tmp(k);%sum(Z_refolded(:,:,cluster_mask,k), 'all');
            % Results(counter).volume_basePeak             = sum(BPC(:));
            Results(counter).basePeak                    = basePeak_mz;%basePeak_mz(k);
            % Results(counter).sample_name                 = matFile_list(k).name;
            Results(counter).sample                      = k;
            Results(counter).sample_feature_ID           = component_number;%counter;
            Results(counter).sample_cluster_ID           = n;
            Results(counter).num_components              = F;

            % Cluster tracking
            cluster_curve_resolution(counter) = component_number;%n + max_cluster;
            counter = counter + 1;
        end
    end

    % Update max cluster id
    max_cluster = max(cluster_curve_resolution);
end

% --- Final trim ---
Results = Results(1:counter-1);
cluster_curve_resolution = cluster_curve_resolution(1:counter-1);

fprintf('Results built: %d total features extracted.\n', numel(Results));
fprintf('Unique base peaks: %d\n', numel(unique([Results.basePeak])));
end