function [merged_borders, Z_tmp, peak_apex, ...
          sum_elution_profile_tmp, groups] = ...
    merge_overlapping_features(borders, data_raw, n)
%MERGE_OVERLAPPING_FEATURES Merge overlapping feature regions.
%
% INPUTS:
%   borders     - cell array:
%                 borders{i}{1} = [x_start x_end]
%                 borders{i}{2} = [y_start y_end]
%   data_raw    - 3D raw data matrix
%   n           - slice/sample index
%
% OUTPUTS:
%   merged_borders          - merged borders for each group
%   Z_tmp                   - merged raw data for each group
%   peak_apex               - [x y] apex position for each group
%   sum_elution_profile_tmp - summed elution profiles
%   groups                  - connected-component group assignment
%
% Overlap between features should be supplied through the
% 'overlap' matrix.

    % ---------------------------------------------------------
    % Find overlapping pairs
    % ---------------------------------------------------------

    [r,c] = find(triu(overlap,1));

    nFeatures = numel(borders);

    % ---------------------------------------------------------
    % Find connected groups
    % ---------------------------------------------------------

    G = graph(r, c, [], nFeatures);
    groups = conncomp(G);

    nGroups = max(groups);

    % ---------------------------------------------------------
    % Extract individual borders
    % ---------------------------------------------------------

    start_1 = cellfun(@(x) x{1}(1), borders);
    end_1   = cellfun(@(x) x{1}(2), borders);

    start_2 = cellfun(@(x) x{2}(1), borders);
    end_2   = cellfun(@(x) x{2}(2), borders);

    % ---------------------------------------------------------
    % Merge borders within each connected group
    % ---------------------------------------------------------

    min_start_1 = accumarray(groups(:), start_1(:), [], @min);
    max_end_1   = accumarray(groups(:), end_1(:),   [], @max);

    min_start_2 = accumarray(groups(:), start_2(:), [], @min);
    max_end_2   = accumarray(groups(:), end_2(:),   [], @max);

    % ---------------------------------------------------------
    % Create output cells
    % ---------------------------------------------------------

    merged_borders = cell(nGroups,1);
    Z_tmp = cell(nGroups,1);
    peak_apex = zeros(nGroups,2);
    sum_elution_profile_tmp = cell(nGroups,1);

    % ---------------------------------------------------------
    % Extract merged regions and find apex
    % ---------------------------------------------------------

    for g = 1:nGroups

        % Merged borders
        border_1 = [min_start_1(g), max_end_1(g)];
        border_2 = [min_start_2(g), max_end_2(g)];

        merged_borders{g} = {
            border_1
            border_2
        };

        % Extract merged raw region
        Z = data_raw( ...
            border_1(1):border_1(2), ...
            border_2(1):border_2(2), ...
            n);

        Z_tmp{g} = Z;

        % -----------------------------------------------------
        % Elution profiles
        % -----------------------------------------------------

        sum_elution_profile_tmp{g} = {
            sum(Z,2)
            sum(Z,1)
        };

        % -----------------------------------------------------
        % Find peak apex
        % -----------------------------------------------------

        [maxIntensity, linearIdx] = max(Z(:));

        [x_local, y_local] = ind2sub(size(Z), linearIdx);

        % Convert local coordinates to original data coordinates
        x_apex = border_1(1) + x_local - 1;
        y_apex = border_2(1) + y_local - 1;

        peak_apex(g,:) = [x_apex, y_apex];

        % fprintf(['Group %d/%d: apex = (%d,%d), ' ...
        %          'intensity = %.4g\n'], ...
        %     g, nGroups, x_apex, y_apex, maxIntensity);

    end

end