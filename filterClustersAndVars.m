function [peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, cosMat, cosMat_combined] = ...
    filterClustersAndVars(peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, cosMat, cosMat_combined, Options)
%FILTERCLUSTERSANDVARS filters clusters{k} and associated variables
%   Keeps only rows that occur more than Options.num_fragments times.

for k = 1:numel(clusters)
    C = clusters{k};
    if isempty(C)
        continue;
    end

    % Find unique rows + counts
    [~, ~, idx] = unique(C, 'rows');
    counts = accumarray(idx, 1);

    % Keep only rows that occur more than Options.num_fragments
    keepIdx = counts(idx) > Options.num_fragments;

    % Apply filtering consistently (check sizes before applying)
    clusters{k}            = C(keepIdx, :);

    if size(peak_borders{k},1) == numel(keepIdx)
        peak_borders{k}    = peak_borders{k}(keepIdx, :);
    end
    if size(sum_elution_profile{k},1) == numel(keepIdx)
        sum_elution_profile{k} = sum_elution_profile{k}(keepIdx, :);
    end
    if size(Z_feature{k},1) == numel(keepIdx)
        Z_feature{k}       = Z_feature{k}(keepIdx, :);
    end
    if size(allPeaks{k},1) == numel(keepIdx)
        allPeaks{k}        = allPeaks{k}(keepIdx, :);
    end
%     if size(cosMat{k},1) == numel(keepIdx)
%         cosMat{k}          = cosMat{k}(keepIdx, keepIdx);
%     end
%     if size(cosMat_combined{k},1) == numel(keepIdx)
%         cosMat_combined{k} = cosMat_combined{k}(keepIdx, keepIdx);
%     end
end
end

