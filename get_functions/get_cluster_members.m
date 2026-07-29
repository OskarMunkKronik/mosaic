function n = get_cluster_members(feature_groups_all, clusters_all_samples, mzTarget, k,Options)
%GET_CLUSTER_MEMBERS  Find all features in the same cluster and sample
%
%   n = GET_CLUSTER_MEMBERS(feature_groups_all, clusters_all_samples, mzTarget, k)
%
%   Inputs:
%       feature_groups_all     – struct array with field 'basePeak' and 'sample'
%       clusters_all_samples   – numeric vector of cluster assignments
%       mzTarget               – target m/z value to match
%       k                      – sample index
%
%   Output:
%       n                      – indices of features in the same cluster and sample
    n = [];
    mz_vec = [feature_groups_all.basePeak];
    sample_vec = [feature_groups_all.sample];
    % Find the closest feature in m/z
    [val, a] = min(abs(mz_vec - mzTarget));
    
    if val/mzTarget * 1e6 < Options.mzTolerance & nnz(feature_groups_all(a).basePeak == mz_vec & sample_vec == k) > 0 
    % Find all features in the same cluster and same sample
    n = find(clusters_all_samples(a) == clusters_all_samples & ...
             sample_vec' == k);
    end 
end
