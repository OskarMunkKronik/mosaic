function run_cluster_chunks(Options, x_range, y_range, mz_range, ...
    Filter, num_clusters, n_samples, chunk_size)



n_digits_samples  = ceil(log10(n_samples ) + 1);
n_digits_clusters = ceil(log10(num_clusters ) + 1);

save_dir = fullfile(Options.dir.results, 'clusters_for_curve_resolution');
q = progressParfor(ceil(num_clusters/chunk_size));

parfor n_cluster = 1:num_clusters
        % n_cluster=vec(n_h);
    send(q,n_cluster)
    % chunk_end = min(chunk_start + chunk_size - 1, num_clusters);
    % cluster_range = chunk_start:chunk_end;

    % fprintf('Processing clusters %d → %d\n', chunk_start, chunk_end);

    %% === STORAGE for this chunk ===
    Z_tmp = cell(1, n_samples);

    %% === LOAD EACH SAMPLE ONCE (PARALLEL) ===
    % q = progressParfor(n_samples);
    for k = 1:n_samples

        % send(q, k);
        sample_str = sprintf(['%0', num2str(n_digits_samples), 'd'], k);
        file_path  = fullfile(save_dir, sprintf('sample_%s.h5', sample_str));


        % Load only once per worker
        Z_local = load_sample_chunk_local( ...
            n_cluster, ...
            file_path, n_digits_clusters);

        % Valid sliced assignment
        Z_tmp(:, k) = Z_local;
    end

    %% === PROCESS EACH CLUSTER (PARALLEL) ===
    % for idx = 1:numel(cluster_range)
    % n_cluster = cluster_range(idx);
    % Z_tmp = Z_chunk(idx, :);
    if all(cellfun(@isempty, Z_tmp))
        continue;
    end
    Z_data_local        = cell(numel(Z_tmp), 2);
    Z_data_local(:, 1)  = Z_tmp(:);
    Z_data_local(:, 2)  = num2cell(1:numel(Z_tmp))';
    Z_tmp_built = build_4D_cluster_cube(Z_data_local, n_cluster, n_samples);
    [W, H, cg, ~, ~] = perform_nnmf_unfold(Z_tmp_built, Options);
    % save_curve_resolution_results( ...
    %     W, H, cg, ...
    %     x_range{n_cluster}, ...
    %     y_range{n_cluster}, ...
    %     mz_range{n_cluster}, ...
    %     n_samples, n_cluster, num_clusters, ...
    %     Options, 'curve_resolution_results');
    save_curve_resolution_results( ...
            [], [], cg,Z_tmp_built,...
            x_range{n_cluster}, ...
            y_range{n_cluster}, ...
            mz_range{n_cluster}, ...
            n_samples, n_cluster, num_clusters, ...
            Options, 'curve_resolution_results');

end

%% === Free memory ===
% clear Z_chunk

end

% end


% helper
% function Z_local = load_sample_chunk_local( ...
%     k, cluster_range, ...
%     save_dir, n_digits_samples, n_digits_clusters, q)
% % Returns ONLY data for this sample (column)
%
%     % Preallocate local column
%     Z_local = cell(numel(cluster_range), 1);
%
%     % Progress (optional)
%     if ~isempty(q)
%         send(q, k);
%     end
%
%     % Build filename
%     sample_str = sprintf(['%0', num2str(n_digits_samples), 'd'], k);
%     file_path  = fullfile(save_dir, sprintf("sample_%s.mat", sample_str));
%
%     if ~exist(file_path, 'file')
%         return;
%     end
%
%     % === FAST: load full sample ONCE ===
%     tmp = load(file_path);
%
%     % === Extract only needed clusters ===
%     for idx = 1:numel(cluster_range)
%
%         n_cluster = cluster_range(idx);
%         cluster_str = sprintf(['cluster_%0', num2str(n_digits_clusters), 'd'], n_cluster);
%
%         if isfield(tmp, cluster_str)
%             s = tmp.(cluster_str);
%             if isfield(s,'Z')
%             Z_local{idx} = s.Z;
%             end
%         end
%     end
%
% end
function Z_local = load_sample_chunk_local( ...
    n_cluster, ...
    file_path, n_digits_clusters)

Z_local = cell(1, 1);
% === Load only requested clusters ===
% for idx = 1:numel(cluster_range)
    % n_cluster   = cluster_range(idx);
    cluster_str = sprintf(['%0', num2str(n_digits_clusters), 'd'], n_cluster);
    group       = sprintf('/cluster_%s', cluster_str);

    % Skip h5info check — just try to read directly, catch empty/missing
    try
        subs        = h5read(file_path, [group '/subs']);
        vals        = h5read(file_path, [group '/vals']);
        tensor_size = h5read(file_path, [group '/tensor_size']);
        Z_local{1} = accumarray(subs, vals, tensor_size);

    catch
        % Empty marker or missing group — leave as []
    end
% end
end