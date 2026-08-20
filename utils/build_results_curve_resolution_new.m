function [Results, cluster_curve_resolution] = ...
    build_results_curve_resolution_new( ...
    W, H, x_range, y_range, mz_range, cg, peak_borders, ...
    mz_ind, mzroi_aug, matFile_list, n_samples, Options)

fprintf('Building Results structure...\n');

% ========================================================================
% General setup
% ========================================================================

num_clusters = numel(mz_ind);

load_dir_curve = fullfile( ...
    Options.dir.results, ...
    'curve_resolution_results');

n_digits_clusters = max(1, ceil(log10(max(num_clusters,1)) + 1));

% ------------------------------------------------------------------------
% Empty result structure
% ------------------------------------------------------------------------

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

% ========================================================================
% FIRST PASS
%
% Load only cg from every existing file.
%
% This lets us determine:
%   1. which clusters exist
%   2. number of components per cluster
%   3. global component offset
%
% This is done before parfor because component_number must be deterministic.
% ========================================================================

cluster_exists = false(num_clusters,1);
n_components = zeros(num_clusters,1);

for n_cluster = 1:num_clusters

    cluster_str = sprintf( ...
        ['%0',num2str(n_digits_clusters),'d'], ...
        n_cluster);

    filename = fullfile( ...
        load_dir_curve, ...
        sprintf('results_%s.mat',cluster_str));

    if ~isfile(filename)
        continue
    end

    S = load(filename,'cg');

    cg_local = S.cg;

    % If 0 means background/unassigned, remove it here.
    component_ids = unique(cg_local(:));
    component_ids = component_ids(component_ids ~= 0);

    cluster_exists(n_cluster) = true;
    n_components(n_cluster) = numel(component_ids);

end

% ========================================================================
% Global component offsets
% ========================================================================

component_offset = zeros(num_clusters,1);

if num_clusters > 1
    component_offset(2:end) = ...
        cumsum(n_components(1:end-1));
end

% Total number of output Results
n_total_results = sum(n_components) * n_samples;

fprintf('Clusters found: %d / %d\n', ...
    sum(cluster_exists), num_clusters);

fprintf('Components: %d\n', sum(n_components));

fprintf('Maximum possible Results entries: %d\n', ...
    n_total_results);

% ========================================================================
% Local output containers
%
% Each parfor iteration produces:
%
%   Results_by_cluster{n_cluster}
%   cluster_ids_by_cluster{n_cluster}
%
% No shared counter is needed.
% ========================================================================

Results_by_cluster = cell(num_clusters,1);
cluster_ids_by_cluster = cell(num_clusters,1);

% ========================================================================
% Parallel processing
% ========================================================================

q = progressParfor(num_clusters);

parfor n_cluster = 1:num_clusters

    send(q,n_cluster);

    % ------------------------------------------------------------
    % Skip missing clusters
    % ------------------------------------------------------------

    if ~cluster_exists(n_cluster)
        Results_by_cluster{n_cluster} = emptyResult([] ,1);
        cluster_ids_by_cluster{n_cluster} = zeros(0,1);
        continue
    end

    % ------------------------------------------------------------
    % Filename
    % ------------------------------------------------------------

    cluster_str = sprintf( ...
        ['%0',num2str(n_digits_clusters),'d'], ...
        n_cluster);

    filename = fullfile( ...
        load_dir_curve, ...
        sprintf('results_%s.mat',cluster_str));

    % ------------------------------------------------------------
    % Load data
    % ------------------------------------------------------------

    S = load(filename,'cg','Z_refolded');

    cg_local = S.cg;
    Z_refolded = S.Z_refolded;

    clear S

    % ------------------------------------------------------------
    % Component IDs
    % ------------------------------------------------------------

    component_ids = unique(cg_local(:));
    component_ids = component_ids(component_ids ~= 0);

    n_components_local = numel(component_ids);

    % ------------------------------------------------------------
    % Peak borders
    % ------------------------------------------------------------

    X = peak_borders{n_cluster}(1,1): ...
        peak_borders{n_cluster}(1,2);

    Y = peak_borders{n_cluster}(2,1): ...
        peak_borders{n_cluster}(2,2);

    nX = numel(X);
    nY = numel(Y);

    % ------------------------------------------------------------
    % m/z indices for this cluster
    % ------------------------------------------------------------

    mz_idx_cluster = mz_ind{n_cluster};

    mz_vals_cluster = mzroi_aug(mz_idx_cluster);

    % ------------------------------------------------------------
    % Number of results for this cluster
    % ------------------------------------------------------------

    nResults_local = n_components_local * n_samples;

    localResults(nResults_local,1) = emptyResult;

    localClusterIDs = zeros(nResults_local,1);

    localCounter = 0;

    % ============================================================
    % COMPONENT LOOP
    % ============================================================

    for component_idx = 1:n_components_local

        component_id = component_ids(component_idx);

        % --------------------------------------------------------
        % Component mask
        % --------------------------------------------------------

        cluster_mask = ...
            cg_local == component_id;

        mz_ind_tmp = ...
            mz_idx_cluster(cluster_mask);

        mz_vals_tmp = ...
            mz_vals_cluster(cluster_mask);

        nMz = numel(mz_ind_tmp);

        if nMz == 0
            continue
        end

        % --------------------------------------------------------
        % Extract component
        %
        % Dimensions:
        %
        %   nX × nY × nMz × n_samples
        % --------------------------------------------------------

        Z_component = ...
            Z_refolded(:,:,cluster_mask,:);

        % --------------------------------------------------------
        % Reshape spatial dimensions
        %
        % P × nMz × n_samples
        %
        % P = nX*nY
        % --------------------------------------------------------

        P = nX * nY;

        Z_2D = reshape( ...
            Z_component, ...
            P, ...
            nMz, ...
            n_samples);

        % ========================================================
        % MASS SPECTRA
        %
        % Sum x/y dimensions.
        %
        % nMz × n_samples
        % ========================================================

        mass_spec = squeeze(sum(Z_2D,1));

        if n_samples == 1
            mass_spec = reshape( ...
                mass_spec, ...
                nMz, ...
                1);
        end

        % ========================================================
        % BASE PEAK m/z
        % ========================================================

        [~,basepeak_idx] = ...
            max(mass_spec,[],1);

        basePeak_mz = ...
            mz_vals_tmp(basepeak_idx);

        basePeak_mz = basePeak_mz(:);

        % ========================================================
        % TOTAL VOLUME
        %
        % Since mass_spec is already summed over x/y:
        %
        % volume = sum over m/z
        % ========================================================

        volume_summed = ...
            sum(mass_spec,1).';

        % ========================================================
        % BASE-PEAK 2D MATRICES
        %
        % Select the base-peak m/z slice for EVERY sample at once.
        %
        % This replaces:
        %
        % for k = 1:n_samples
        %     BPC = Z_component(:,:,bp_idx,k);
        % end
        %
        % ========================================================

        % Reshape to:
        %
        % P × (nMz*n_samples)

        Z_flat = reshape( ...
            Z_2D, ...
            P, ...
            nMz*n_samples);

        % Columns containing the selected base peak
        %
        % For sample k:
        %
        % column = basepeak_idx(k) + (k-1)*nMz

        sample_offset = ...
            (0:n_samples-1) * nMz;

        basepeak_columns = ...
            basepeak_idx + sample_offset;

        % Extract all BPCs at once:
        %
        % P × n_samples

        BPC_all = ...
            Z_flat(:,basepeak_columns);

        % ========================================================
        % Find maximum spatial position for every sample
        % ========================================================

        [~,linear_idx] = ...
            max(BPC_all,[],1);

        linear_idx = linear_idx(:);

        [r,c] = ind2sub( ...
            [nX,nY], ...
            linear_idx);

        r = r(:);
        c = c(:);

        % Convert to actual X/Y values
        Rt_x = X(r);
        Rt_y = Y(c);

        % ========================================================
        % GLOBAL COMPONENT ID
        % ========================================================

        component_number = ...
            component_offset(n_cluster) + component_id;

        % --------------------------------------------------------
        % Save component
        % --------------------------------------------------------

        save_curve_resolution_results( ...
            [],[], ...
            component_id, ...
            Z_component, ...
            X,Y, ...
            mz_ind_tmp, ...
            n_samples, ...
            component_number, ...
            num_clusters, ...
            Options, ...
            'component_cluster');

        % ========================================================
        % BUILD SPARSE MASS SPECTRA
        %
        % Instead of constructing a sparse vector inside the
        % sample loop, build ALL samples simultaneously.
        %
        % Dimensions:
        %
        % length(mzroi_aug) × n_samples
        % ========================================================

        sparse_rows = repmat( ...
            mz_ind_tmp(:), ...
            n_samples,1);

        sparse_cols = kron( ...
            (1:n_samples).', ...
            ones(nMz,1));

        sparse_values = mass_spec(:);

        massSpec_all = sparse( ...
            sparse_rows, ...
            sparse_cols, ...
            sparse_values, ...
            numel(mzroi_aug), ...
            n_samples);

        % ========================================================
        % RESULT INDICES
        % ========================================================

        result_idx = ...
            localCounter + (1:n_samples);

        % ========================================================
        % Fill structure array
        %
        % These assignments are vectorized over samples.
        % ========================================================

        localResults(result_idx).Rt = ...
            num2cell([Rt_x.'; Rt_y.'],1).';

        localResults(result_idx).elution1 = ...
            repmat({X},n_samples,1);

        localResults(result_idx).elution2 = ...
            repmat({Y},n_samples,1);

        localResults(result_idx).massSpec = ...
            num2cell(massSpec_all,1).';

        localResults(result_idx).volume_summed = ...
            num2cell(volume_summed);

        localResults(result_idx).basePeak = ...
            num2cell(basePeak_mz);

        localResults(result_idx).sample = ...
            num2cell((1:n_samples).');

        localResults(result_idx).sample_feature_ID = ...
            repmat({component_number},n_samples,1);

        localResults(result_idx).sample_cluster_ID = ...
            repmat({component_id},n_samples,1);

        localResults(result_idx).num_components = ...
            repmat({n_components_local},n_samples,1);

        localClusterIDs(result_idx) = ...
            component_number;

        % --------------------------------------------------------
        % Advance local counter
        % --------------------------------------------------------

        localCounter = ...
            localCounter + n_samples;

        % --------------------------------------------------------
        % Free large temporary arrays
        % --------------------------------------------------------

        % clear Z_component Z_2D Z_flat BPC_all
        % clear mass_spec massSpec_all

    end

    % ============================================================
    % Trim local arrays
    % ============================================================

    if localCounter < nResults_local

        localResults = ...
            localResults(1:localCounter);

        localClusterIDs = ...
            localClusterIDs(1:localCounter);

    end

    Results_by_cluster{n_cluster} = ...
        localResults;

    cluster_ids_by_cluster{n_cluster} = ...
        localClusterIDs;

end

% ========================================================================
% Concatenate results
% ========================================================================

fprintf('Combining cluster results...\n');

valid_clusters = ...
    find(cluster_exists);

if isempty(valid_clusters)

    Results = emptyResult([]);
    cluster_curve_resolution = zeros(0,1);

    fprintf('No results found.\n');
    return

end

% Concatenate struct arrays
Results = vertcat( ...
    Results_by_cluster{valid_clusters});

cluster_curve_resolution = vertcat( ...
    cluster_ids_by_cluster{valid_clusters});

% ========================================================================
% Final consistency check
% ========================================================================

if numel(Results) ~= numel(cluster_curve_resolution)

    error([ ...
        'Internal indexing error: Results has %d entries, ' ...
        'cluster_curve_resolution has %d entries.'], ...
        numel(Results), ...
        numel(cluster_curve_resolution));

end

fprintf( ...
    'Results built: %d total features extracted.\n', ...
    numel(Results));

if ~isempty(Results)

    fprintf( ...
        'Unique base peaks: %d\n', ...
        numel(unique([Results.basePeak])));

end

end