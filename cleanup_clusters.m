function [peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, cosMat, cosMat_combined] = ...
    cleanup_clusters(k, peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, cosMat, cosMat_combined)
% CLEANUP_CLUSTERS Remove empty/zeroed entries after cluster filtering
%
% Usage:
%   [peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, cosMat, cosMat_combined] = ...
%       cleanup_clusters(k, peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, cosMat, cosMat_combined)
%
% Input:
%   k                - dataset index
%   peak_borders     - cell array of peak borders
%   sum_elution_profile - cell array of elution profiles
%   Z_feature        - cell array of features
%   allPeaks         - numeric array of peaks
%   clusters         - cluster assignments
%   cosMat           - cell array of cosine similarity matrices {1}, {2}
%   cosMat_combined  - combined cosine similarity matrix
%
% Output:
%   Same variables with entries removed where clusters{k} == 0
%
% Author: Oskar Munk Kronik, oskarmunkkronik@gmail.com
% Date: 17-09-2025
if size(allPeaks,2) == 5
    allPeaks{k} = cat(2,allPeaks{k},zeros(length(clusters{k}),1));
else
    allPeaks{k}(:,6) = 0;
end
for n_cluster = 1:length(clusters{k})
    ind_cluster = find(clusters{k} == n_cluster);
    [~,ind_max]= max(allPeaks{k}(ind_cluster,5));
    % allPeaks{k}(ind_cluster,[2,3]) =  repmat(allPeaks{k}(ind_cluster(ind_max),[2,3]),[length(ind_cluster),1]);
    allPeaks{k}(ind_cluster(ind_max),6) = n_cluster;
    if any(groupcounts(allPeaks{k}(ind_cluster,4)) > 1)

        [~,~,u_vec] = unique(allPeaks{k}(ind_cluster,4));
        max_val = accumarray(u_vec,allPeaks{k}(ind_cluster,5),[],@max);
        %         [~, max_ind] = max(allPeaks{k}(ind_cluster,5));
        % fill in max volume Rt for all in that feature group

        %         ind_rt = ind_cluster(ismember(allPeaks{k}(ind_cluster,5),max_val));
        ind_cluster = ind_cluster(~ismember(allPeaks{k}(ind_cluster,5),max_val));

        % Mark for deletion
        peak_borders{k}(ind_cluster) = {[]};
        sum_elution_profile{k}(ind_cluster) = {[]};
        Z_feature{k}(ind_cluster) = {[]};
        allPeaks{k}(ind_cluster,:) = 0;
        clusters{k}(ind_cluster) = 0;

        %         for dim = 1:2
        % %             cosMat{k}(ind_cluster,ind_cluster,dim) = 0;
        %         end
        %         cosMat_combined{k}(ind_cluster,ind_cluster) = 0;
    end
end
% indices marked for deletion
toRemove = clusters{k} == 0;

% --- cell arrays ---
peak_borders{k}(toRemove)        = [];
sum_elution_profile{k}(toRemove) = [];
Z_feature{k}(toRemove)           = [];

% --- numeric arrays ---
allPeaks{k}(toRemove,:) = [];
clusters{k}(toRemove)   = [];

% --- similarity matrices ---
%     cosMat{k}(toRemove,:,:) = [];
%     cosMat{k}(:,toRemove,:) = [];
%     cosMat{k}(toRemove,:,2) = [];
%     cosMat{k}(:,toRemove,2) = [];

%     cosMat_combined{k}(toRemove,:) = [];
%     cosMat_combined{k}(:,toRemove) = [];
end
