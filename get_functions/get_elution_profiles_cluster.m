function Z_aligned = get_elution_profiles_cluster(Z_feature, peak_borders, clusters, points, k)
%ALIGN_PEAK_CLUSTERS Aligns 2D elution maps across features in each cluster
%
%   Z_aligned = align_peak_clusters(Z_feature, peak_borders, clusters, points, k, mzroi_aug, target_mz)
%
%   Inputs:
%       Z_feature     - cell array of 2D elution maps (e.g. Z_feature{k}{i})
%       peak_borders  - nested cell array defining peak borders per feature
%       clusters      - cell array of cluster indices per dataset
%       points        - struct array with field .Rt (retention times)
%       k             - dataset index (integer)
%       mzroi_aug     - vector of m/z ROI centers
%       target_mz     - scalar m/z target to locate cluster (optional)
%
%   Output:
%       Z_aligned     - cell array of aligned 3D cubes (each cluster's Z_tmp)
%
%   Notes:
%       Pads all peaks in each cluster to common min/max retention borders.
%       Keeps the original elution intensity values in correct relative positions.

% -------------------------------------------------------------
% Step 1: Prepare retention time and find non-empty clusters
% -------------------------------------------------------------
Rt_mat = [points.Rt]';
nDim = 2;

flag = zeros(max(clusters{k}), size(Rt_mat, 2));
for dim = 1:nDim
    flag(:, dim) = accumarray(clusters{k}, Rt_mat(:, dim), [], @max) - ...
                   accumarray(clusters{k}, Rt_mat(:, dim), [], @min);
end
flag = find(sum(flag, 2) > 0);

% % Optionally locate target m/z for debugging or filtering
% if nargin >= 7 && ~isempty(target_mz)
%     [~, a] = min(abs(mzroi_aug - target_mz));
%     fprintf('Closest m/z index to %.4f found at index %d\n', target_mz, a);
% end

% -------------------------------------------------------------
% Step 2: Build aligned cubes for each cluster
% -------------------------------------------------------------
Z_aligned = cell(length(peak_borders{k}), 1);

for n = 1:length(flag)
    vec_tmp = find(clusters{k} == flag(n));

    % --- Collect peak borders ---
    border_1D = zeros(length(vec_tmp), nDim);
    border_2D = zeros(length(vec_tmp), nDim);
    for nn = 1:length(vec_tmp)
        border_1D(nn, :) = peak_borders{k}{vec_tmp(nn)}{1};
        border_2D(nn, :) = peak_borders{k}{vec_tmp(nn)}{2};
    end

    % --- Global border extents ---
    min_border_1D = min(border_1D, [], 1);
    max_border_1D = max(border_1D, [], 1);
    min_border_2D = min(border_2D, [], 1);
    max_border_2D = max(border_2D, [], 1);

    % --- Initialize padded elution cube ---
    Z_tmp = zeros( ...
        max_border_1D(2), ...
        max_border_2D(2), ...
        length(vec_tmp) );

    % --- Fill each feature's data into aligned cube ---
    for nn = 1:length(vec_tmp)
        Z_nn = Z_feature{k}{vec_tmp(nn)};

        range_1D = border_1D(nn,1):border_1D(nn,2);
        range_2D = border_2D(nn,1):border_2D(nn,2);

        Z_tmp(range_1D, range_2D, nn) = Z_nn;
    end

    % --- Crop to cluster boundaries ---
    Z_tmp = Z_tmp(min_border_1D(1):max_border_1D(2), ...
                  min_border_2D(1):max_border_2D(2), :);

    Z_aligned{flag(n)} = Z_tmp;
end
end
