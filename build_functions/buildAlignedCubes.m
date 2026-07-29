function Z_aligned = buildAlignedCubes(peak_borders, clusters, Z_feature, flag, nDim, k)
% buildAlignedCubes
% -------------------------------------------------------------
% Build aligned elution cubes for each cluster in sample k.
%
% INPUTS:
%   peak_borders : cell array {k}{feature}{dim}, containing
%                  start and end indices for each feature in both dimensions.
%   clusters     : cell array {k}, containing cluster labels for each feature.
%   Z_feature    : cell array {k}{feature}, containing each feature's 2D data.
%   flag         : vector of unique cluster identifiers to process.
%   nDim         : number of chromatographic dimensions (usually 2).
%   k            : index of the current sample or dataset.
%
% OUTPUT:
%   Z_aligned    : cell array where each entry corresponds to an aligned
%                  3D cube (1st and 2nd dimensions = elution, 3rd = features).
%
% EXAMPLE:
%   Z_aligned = buildAlignedCubes(peak_borders, clusters, Z_feature, unique(clusters{k}), 2, k);
%
% -------------------------------------------------------------

% Preallocate output
Z_aligned = cell(length(peak_borders{k}), 1);

for n = 1:length(flag)
    % --- Identify all features belonging to this cluster
    vec_tmp = find(clusters{k} == flag(n));
    if isempty(vec_tmp), continue; end

    % --- Collect 1D and 2D peak borders for all features
    border_1D = zeros(length(vec_tmp), nDim);
    border_2D = zeros(length(vec_tmp), nDim);
    for nn = 1:length(vec_tmp)
        border_1D(nn, :) = peak_borders{k}{vec_tmp(nn)}{1};
        border_2D(nn, :) = peak_borders{k}{vec_tmp(nn)}{2};
    end

    % --- Determine global cluster border extents
    min_border_1D = min(border_1D, [], 1);
    max_border_1D = max(border_1D, [], 1);
    min_border_2D = min(border_2D, [], 1);
    max_border_2D = max(border_2D, [], 1);

    % --- Initialize padded elution cube
    Z_tmp = zeros(max_border_1D(2), max_border_2D(2), length(vec_tmp));

    % --- Fill each feature's data into the aligned cube
    for nn = 1:length(vec_tmp)
        Z_nn = Z_feature{k}{vec_tmp(nn)};

        range_1D = border_1D(nn,1):border_1D(nn,2);
        range_2D = border_2D(nn,1):border_2D(nn,2);

        Z_tmp(range_1D, range_2D, nn) = Z_nn;
    end

    % --- Crop cube to the cluster boundaries
    Z_tmp = Z_tmp(min_border_1D(1):max_border_1D(2), ...
                  min_border_2D(1):max_border_2D(2), :);

    % --- Store in output
    Z_aligned{flag(n)} = Z_tmp;
end
end
