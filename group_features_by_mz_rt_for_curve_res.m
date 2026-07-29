function [clusters_final, mz_merged, peak_borders_final] = group_features_by_mz_rt_for_curve_res(feature_groups_all, clusters_all_samples, Options)
%GROUP_FEATURES_BY_MZ_RT Groups features first by RT similarity and then by overlapping m/z
%
% Inputs:
%   feature_groups_all  : array of structs, each with fields:
%                         - Rt          : [2 x 1] retention times
%                         - massSpec    : vector of intensities
%                         - elution1/2  : elution profiles (for 1D/2D)
%   clusters_all_samples: numeric cluster assignment vector
%   Options             : struct with field .Rtdev (RT deviation threshold)
%
% Outputs:
%   clusters_final      : final cluster assignments after merging
%   mz_merged           : merged m/z index sets for each cluster
%   peak_borders_final  : merged 1D/2D peak borders for each cluster

%% ------------------------------------------------------------------------
% --- 1. Collect RT and elution borders ----------------------------------
nClustersAll = max(clusters_all_samples);
Rt_mat = cell2mat(arrayfun(@(x) x.Rt, feature_groups_all, 'UniformOutput', false));

% Median RT per cluster
rt_mat_cluster = zeros(2, nClustersAll);
for dim = 1:2
    rt_mat_cluster(dim,:) = round(accumarray(clusters_all_samples, Rt_mat(dim,:), [nClustersAll,1], @median, NaN));
end

% Collect elution borders
peak_borders_tmp = cell(2,1);
peak_borders_tmp{1}(:,1) = arrayfun(@(x) find(x.elution1 > 0,1,'first'), feature_groups_all, 'UniformOutput', false);
peak_borders_tmp{1}(:,2) = arrayfun(@(x) find(x.elution1 > 0,1,'last'),  feature_groups_all, 'UniformOutput', false);
peak_borders_tmp{2}(:,1) = arrayfun(@(x) find(x.elution2 > 0,1,'first'), feature_groups_all, 'UniformOutput', false);
peak_borders_tmp{2}(:,2) = arrayfun(@(x) find(x.elution2 > 0,1,'last'),  feature_groups_all, 'UniformOutput', false);

%% ------------------------------------------------------------------------
% --- 2. Group clusters by RT proximity ----------------------------------
idx = 1:nClustersAll;
refs = idx(1);
i_vec = idx(2:end);
max_dist = Options.Rtdev(:);

clusters_rt = ones(nClustersAll,1);
clusters_old = clusters_all_samples(:);
u_clusters = unique(clusters_all_samples);

for i = i_vec
    validRefs = all(abs(rt_mat_cluster(:,i) - rt_mat_cluster(:,refs)) <= max_dist, 1);
    passIdx = find(validRefs, 1, 'first');
    if any(validRefs)
        clusters_rt(i) = clusters_rt(refs(passIdx));
        clusters_old(i) = u_clusters(refs(passIdx));
    else
        refs(end+1) = i;
        clusters_rt(i) = numel(refs);
        clusters_old(i) = u_clusters(i);
    end
end

%% ------------------------------------------------------------------------
% --- 3. Collect massSpec indices per RT cluster -------------------------
nRtClusters = max(clusters_rt);
[massSpecVals, mz_ind, peak_borders] = deal(cell(nRtClusters,1));

for n_cluster = 1:nRtClusters
    ind_cluster = any(clusters_all_samples == clusters_old(clusters_rt == n_cluster)', 2);

    active_mass_idx = arrayfun(@(x) find(x.massSpec > 0), feature_groups_all(ind_cluster), 'UniformOutput', false);
    massSpecVals{n_cluster} = active_mass_idx;
    mz_ind{n_cluster} = active_mass_idx; % not merged yet

    for dim = 1:2
        peak_borders{n_cluster}{dim} = cell2mat(peak_borders_tmp{dim}(ind_cluster,:))';
    end
end

%% ------------------------------------------------------------------------
% --- 4. Merge RT clusters sharing m/z indices ---------------------------
[peak_borders_final, clusters_final, mz_merged] = deal(cell(nRtClusters,1));
max_cluster = 0;

for n_cluster = 1:nRtClusters
    nSub = numel(mz_ind{n_cluster});
    if nSub < 2
        clusters_final{n_cluster} = max_cluster + 1;
        tmp = zeros(2);
        for dim = 1:2
            tmp(dim,:) = [peak_borders{n_cluster}{dim}(1), peak_borders{n_cluster}{dim}(2)];
        end
        peak_borders_final{n_cluster,1} = tmp;
        % --- Merge m/z indices per component ---
         mz_merged{n_cluster} =mz_ind{n_cluster}{:};
        max_cluster = max_cluster + 1;
        
        continue;
    end

    % --- Build adjacency matrix for overlapping m/z ---
    match = false(nSub);
    for j = 1:nSub
        for k = j+1:nSub
            if any(ismember(mz_ind{n_cluster}{j}, mz_ind{n_cluster}{k}))
                match(j,k) = true;
                match(k,j) = true;
            end
        end
    end

    % --- Build graph & find connected components ---
    G = graph(match);
    clusters_tmp = conncomp(G)';

    clusters_final{n_cluster} = clusters_tmp + max_cluster;

    % --- Merge m/z indices per component ---
    % mz_merged{n_cluster} = arrayfun(@(c) unique(vertcat(mz_ind{n_cluster}{clusters_tmp == c})), ...
    %     1:max(clusters_tmp), 'UniformOutput', false)';
    mz_merged{n_cluster} = arrayfun(@(c) unique(vertcat(mz_ind{n_cluster}{clusters_tmp == c})), ...
    1:max(clusters_tmp), 'UniformOutput', false)';


    % --- Merge peak borders (min start / max end per dim) ---
    for c = 1:max(clusters_tmp)
        tmp = zeros(2);
        for dim = 1:2
            starts = peak_borders{n_cluster}{dim}(1,clusters_tmp == c);
            ends   = peak_borders{n_cluster}{dim}(2,clusters_tmp == c);
            tmp(dim,:) = [floor(min(starts)), ceil(max(ends))];
        end
        peak_borders_final{n_cluster}{c,1} = tmp;
    end

    max_cluster = max(clusters_final{n_cluster});
end

% Flatten cell-of-cells
% clusters_final = vertcat(clusters_final{:});
% mz_merged = vertcat(mz_merged{:});
% peak_borders_final = vertcat(peak_borders_final{:});

% fprintf('RT grouping reduced %d clusters → %d RT groups → %d final m/z groups.\n', ...
%     nClustersAll, nRtClusters, numel(clusters_final));

end
