
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