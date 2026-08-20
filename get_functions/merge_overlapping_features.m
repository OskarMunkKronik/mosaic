function [merged_borders, Z_tmp, peak_apex, ...
    sum_elution_profile_tmp, max_intensity] = ...
    merge_overlapping_features(overlap, borders, data_raw, n)

% -------------------------------------------------------------
% Number of features
% -------------------------------------------------------------

nFeatures = numel(borders);

% -------------------------------------------------------------
% Find overlapping pairs
% -------------------------------------------------------------

[r,c] = find(triu(overlap,1));

% -------------------------------------------------------------
% Connected components
%
% Every non-overlapping feature is also its own group.
% -------------------------------------------------------------

G = graph(r,c,[],nFeatures);
groups = conncomp(G);

nGroups = max(groups);

% -------------------------------------------------------------
% Size of original data
% -------------------------------------------------------------

sz = size(data_raw);

if numel(sz) < 3
    sz(3) = 1;
end

% -------------------------------------------------------------
% Original borders
% -------------------------------------------------------------

start_1 = cellfun(@(x) x{1}(1), borders);
end_1   = cellfun(@(x) x{1}(2), borders);

start_2 = cellfun(@(x) x{2}(1), borders);
end_2   = cellfun(@(x) x{2}(2), borders);

% -------------------------------------------------------------
% Merged borders
% -------------------------------------------------------------

min_start_1 = accumarray( ...
    groups(:), start_1(:), [], @min);

max_end_1 = accumarray( ...
    groups(:), end_1(:), [], @max);

min_start_2 = accumarray( ...
    groups(:), start_2(:), [], @min);

max_end_2 = accumarray( ...
    groups(:), end_2(:), [], @max);

% -------------------------------------------------------------
% Allocate
% -------------------------------------------------------------

merged_borders = cell(nGroups,1);
Z_tmp = cell(nGroups,1);

sum_elution_profile_tmp = cell(nGroups,1);

peak_apex = zeros(nGroups,2);
max_intensity = zeros(nGroups,1);

% =============================================================
% Process groups
% =============================================================

for g = 1:nGroups

    % ---------------------------------------------------------
    % Merged borders
    % ---------------------------------------------------------

    border_1 = [
        min_start_1(g), ...
        max_end_1(g)
    ];

    border_2 = [
        min_start_2(g), ...
        max_end_2(g)
    ];

    merged_borders{g} = {
        border_1
        border_2
    };

    % ---------------------------------------------------------
    % Extract merged region
    % ---------------------------------------------------------

    Z = data_raw( ...
        border_1(1):border_1(2), ...
        border_2(1):border_2(2), ...
        n);

    Z_tmp{g} = Z;

    % ---------------------------------------------------------
    % Full-size elution profiles
    %
    % Same format as in peak_detection_chunck:
    %
    %   {1} = dimension 1, column vector of size sz(1)
    %   {2} = dimension 2, column vector of size sz(2)
    % ---------------------------------------------------------

    profile_1 = zeros(sz(1),1);
    profile_2 = zeros(sz(2),1);

    range_1 = border_1(1):border_1(2);
    range_2 = border_2(1):border_2(2);

    profile_1(range_1) = sum(Z,2);

    % sum(Z,1) is 1 x N, so transpose it
    profile_2(range_2) = sum(Z,1).';

    sum_elution_profile_tmp{g} = {
        profile_1
        profile_2
    };

    % ---------------------------------------------------------
    % Maximum intensity
    % ---------------------------------------------------------

    [max_intensity(g),linear_idx] = max(Z(:));

    [x_local,y_local] = ind2sub(size(Z),linear_idx);

    % ---------------------------------------------------------
    % Convert local apex coordinates to global coordinates
    % ---------------------------------------------------------

    peak_apex(g,:) = [
        border_1(1) + x_local - 1, ...
        border_2(1) + y_local - 1
    ];

end

end