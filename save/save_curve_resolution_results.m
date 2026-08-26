function save_curve_resolution_results(W,H,cg,Z_refolded, x_range, y_range, mz_range, ...
                              n_samples, n_cluster, num_clusters, ...
                              Options,suffix)
%SAVE_CLUSTER_RESULTS Save NNMF cluster results with structured folder names.
%
%   save_cluster_results(W, H, x_range, y_range, mz_range, ...
%                        n_samples, n_cluster, num_clusters, Options, sample_str)
%
%   Creates zero-padded cluster folders and saves NNMF results per cluster
%   and sample.
%
%   Inputs:
%       W, H            - Factor matrices for the current cluster
%       x_range, y_range, mz_range - Range indices for this cluster
%       n_samples       - Number of samples
%       n_cluster       - Current cluster index
%       num_clusters    - Total number of clusters (for padding)
%       Options         - Struct containing Options.dir.results
%       sample_str      - String or ID for the sample (used in filename)
%
%   Example:
%       save_cluster_results(W, H, x_range, y_range, mz_range, 5, 3, 25, Options, 'SampleA')
    
    % === Setup ===
    n_digits_clusters = ceil(log10(num_clusters ) + 1);
    save_dir = fullfile(Options.Paths.save2mat, suffix);
   

    % Zero-padded cluster folder name
    cluster_str = sprintf(['%0', num2str(n_digits_clusters), 'd'], n_cluster);
    % save_dir_cluster = fullfile(save_dir, cluster_str);

    % Ensure directory exists
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    % Output filename (zero-padded and structured)
    filename = sprintf("results_%s.mat", cluster_str);

    % === Save ==
    save(fullfile(save_dir, filename), ...
         "W","H","cg", "x_range", "y_range", "mz_range", ...
         "n_samples", "n_cluster","Z_refolded",'-nocompression');

    % fprintf('✅ Saved cluster %s results to: %s\n', cluster_str, save_dir_cluster);
end
