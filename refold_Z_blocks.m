function Z_refolded = refold_Z_blocks(W, H, x_range, y_range, mz_range, n_samples, n_cluster)
%REFOLD_Z_BLOCKS Reconstruct and refold Z tensors from factorized components W and H.
%
%   Z_refolded = REFOLD_Z_BLOCKS(W, H, x_range, y_range, mz_range, n_samples)
%
%   Inputs:
%       W           - Left factor matrix from decomposition (e.g., from NMF)
%       H           - Right factor matrix (matching columns of W)
%       x_range     - Cell array of X-dimension index ranges for each cluster
%       y_range     - Cell array of Y-dimension index ranges for each cluster
%       mz_range    - Cell array of m/z index ranges for each cluster
%       n_samples   - Number of samples (fourth dimension)
%
%   Output:
%       Z_refolded  - Cell array where each cell {n_cluster} is a 4D cube:
%                       (x, y, m/z, sample)
%
%   Notes:
%       • Automatically computes sizes from x/y/mz ranges.
%       • Uses parallelization (parfor) for speed.
%       • Each cube is built slice-by-slice and permuted correctly.

    % === Compute per-cluster dimensions ===
    % num_clusters = numel(x_range);
    d1 = numel(x_range);
    d2 = numel(y_range);
    d3 = numel(mz_range);
   
    Z_refolded = permute( ...
        reshape(W * H, [d2, d1, n_samples, d3]), ...
        [2, 1, 4, 3]);

end
