function [mz_ind,peak_borders,gc,grps,Rt,active] = get_mz_ind_for_curve_resolution(d,mzTarget,clusters_all_samples,feature_groups_all,Options,mass_spectrum,active)
% Find all signals within ppm window
% [~,a] = min(abs(mzroi_aug-mzTarget));
% ind = d(:,1) == mzroi_aug(a);
% ppmError = abs(d(:,1) - mzTarget) ./ mzTarget * 1e6;
% ind = ppmError < Options.ppm_dev;
% a = unique(vertcat(mass_spectrum{ind}));
nFeatures = length(mass_spectrum);
ind = false(nFeatures,1);
ppm_dev = Options.ppm_dev;
vec_features = find(active);
% ind = false(nFeatures,1);
% ind = vec_features;
%%
ind(vec_features) = cellfun(@(mz_spec) ...
    any(any(abs(mz_spec - mzTarget) ./ mzTarget * 1e6 < ppm_dev, 2)), ...
    mass_spectrum(vec_features));
%%
% ind = cell2mat(cellfun(@(x) ...
%     any(x == a,'all'), ...
%     mass_spectrum,'UniformOutput',false))>0;


[peak_borders,mz_ind,Rt] = deal(cell(1));
% Find dominant cluster among those matches
[gc, grps] = groupcounts(clusters_all_samples(ind));
if ~any(ind)
    return
end 
active(ind) = false;
% [~, max_group] = max(gc);

% choose_cluster = grps;

% Get all samples in that cluster
[r,c] = find(clusters_all_samples == grps');
% ind_cluster = any(clusters_all_samples == grps',2);

%%
% Rt_mat = cell2mat(arrayfun(@(x) x.Rt, feature_groups_all(r), 'UniformOutput', false));
Rt_mat = [feature_groups_all(r).Rt];
[rt_mat_cluster_min,rt_mat_cluster_max,rt_mat_cluster] = deal(zeros([2,length(grps)]));
for dim = 1:2
    rt_mat_cluster(dim,:) = accumarray(c,Rt_mat(dim,:),[],@median);
    rt_mat_cluster_min(dim,:) = accumarray(c,Rt_mat(dim,:),[],@min);
    rt_mat_cluster_max(dim,:) = accumarray(c,Rt_mat(dim,:),[],@max);

end


idx = 1:length(grps);
refs = idx(1);
i_vec = idx(2:end);%2:n;
max_dist = Options.Rtdev';
clusters = ones(length(grps),1);
for i = i_vec
    assigned = false;
    % Differences to all references at once

    % Keep only refs where both coordinates are within max_dist
    validRefs =all(abs(rt_mat_cluster(:,i) - rt_mat_cluster(:,refs)) <= max_dist, 1);
    passIdx = find(validRefs,1,'first');
    if any(validRefs)
        clusters(i) = passIdx; % or index in refs depending on logic
        assigned = true;
    end

    if ~assigned
        refs(end+1) = i; % new cluster reference
        clusters(i) = length(refs);
        i_vec = i_vec(~(ismember(i_vec,1:length(refs))));
    end
end


for n_cluster = 1:max(clusters)

    ind_cluster_2 = any(clusters_all_samples == grps(clusters == n_cluster)',2);

    % Collect unique active MS features
    massSpecVals{n_cluster} = arrayfun(@(x) find(x.massSpec > 0), feature_groups_all(ind_cluster_2), 'UniformOutput', false);

    %1D
    peak_borders{n_cluster}(1,1) =  floor(median(cell2mat(arrayfun(@(x) find(x.elution1 > 0,1,'first'), feature_groups_all(ind_cluster_2), 'UniformOutput', false))));
    peak_borders{n_cluster}(1,2) =  ceil(median(cell2mat(arrayfun(@(x) find(x.elution1 > 0,1,'last'), feature_groups_all(ind_cluster_2), 'UniformOutput', false))));

    %2D
    peak_borders{n_cluster}(2,1) =  floor(median(cell2mat(arrayfun(@(x) find(x.elution2 > 0,1,'first'), feature_groups_all(ind_cluster_2), 'UniformOutput', false))));
    peak_borders{n_cluster}(2,2) =  ceil(median(cell2mat(arrayfun(@(x) find(x.elution2 > 0,1,'last'), feature_groups_all(ind_cluster_2), 'UniformOutput', false))));

    Rt{n_cluster} = round(median(Rt_mat(:,c == n_cluster),2));
    mz_ind{n_cluster} = unique(vertcat(massSpecVals{n_cluster}{:}));
end


end