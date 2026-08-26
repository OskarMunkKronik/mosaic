% function [feature_table, feature_groups_all, clusters_all_samples,clusters_all_samples_before_duplicate_removal] = cleanup_clusters_feature_table(feature_groups_all,clusters_all_samples,feature_table, Options)
% 
%     cos_thresh = 1-Options.Clustering.cutoff;
%     % Assign same Rt to all with the same cluster: 
%     feature_table_old = feature_table;
%     %%
%     groups = feature_table(:,2);
%     values = feature_table(:,6);
% 
%     maxIdx = accumarray(groups, (1:size(feature_table,1))', [], ...
%         @(idx) idx(find(feature_table(idx,6) == max(feature_table(idx,6)), 1)));
%     mz_median = feature_table(maxIdx,1);
%     Rt_median(:,1:2) = feature_table(maxIdx,[3,4]);
%     %%
%     % feature_table(any(feature_table(:,[1,3,4]) == 0,2),[1,3,4]) =nan;
%     % Rt_median(:,1) = accumarray(feature_table(:,2),feature_table(:,3),[],@(x) median(x,'omitnan'));
%     % Rt_median(:,2) = accumarray(feature_table(:,2),feature_table(:,3),[],@(x) median(x,'omitnan'));
%     Rt_median = Rt_median(feature_table(:,2),:);
%     % 
%     % mz_median = accumarray(feature_table(:,2),feature_table(:,1),[],@(x) median(x,'omitnan'));
%     mz_median = mz_median(feature_table(:,2));
% 
%     % Sort table by mz, RT, sample
%     feature_table(:,1) =  mz_median ;
%     feature_table(:,3:4) = Rt_median(:,:);
%     % [feature_table, sort_idx] = sortrows(feature_table, [1,3,4], "ascend");
% 
%     [feature_table, sort_idx] = sortrows(feature_table, [1,3,4], "ascend");
%     clusters_all_samples_before_duplicate_removal = clusters_all_samples(sort_idx);
%     feature_table_old = feature_table_old(sort_idx,:);
%     if ~isempty(feature_groups_all)
%     feature_groups_all = feature_groups_all(sort_idx);
%     end
% 
%     % Rt_median = Rt_median(sort_idx,:);
%     % mz_median = mz_median(sort_idx);
% 
%     N = size(feature_table,1);
%     clusters = zeros(1,N);
%     clusters(1) = 1;
% 
% 
%     % -------------------------------------------------------------
%     % Build clusters based on RT proximity
%     % -------------------------------------------------------------
% 
%     m = 1;
%     wU = [(1 + Options.ppm_dev*1e-6) * feature_table(m,1), Options.Rtdev];
%     % wU = [0, Options.Clustering.distance_componentization];
%     for n = 2:N
%         %   % --- Compute cosine similarity -
%         % ref = full(feature_groups_all(n).massSpec);
%         % cand = full(feature_groups_all(m).massSpec);
%         % cosine_sim = dot(cand , ref) / (norm(cand) * norm(ref) + eps);
% 
%         if all(abs(feature_table(n,[1,3,4]) - feature_table(m,[1,3,4])) <= wU) %all( feature_table(n,3:4) - feature_table(m,3:4) <= wU(2:3) & cosine_sim > cos_thresh)% 
%        % all(feature_table(n,1) - feature_table(m,1) == wU(1) & feature_table(n,3:4) - feature_table(m,3:4) <= wU(2:3) & cosine_sim > cos_thresh)
%             %    clf 
%         %     stem(ref,'Marker','none')
%         %     hold on
%         % stem(-cand,'Marker','none')
% 
%             continue
%         end
%         m = n;
%         % max_dist_rt = Rt_median(m,:) + Options.Clustering.distance_componentization;
%         wU(1) = (1 + Options.ppm_dev*1e-6) * feature_table(m,1);
%         clusters(1,n) = 1;
%     end
% 
%     clusters = cumsum(clusters);
% 
%     % -------------------------------------------------------------
%     % PROCESS clusters without deleting inside the loop
%     % -------------------------------------------------------------
%     max_clust = max(clusters);
% 
%     q = progressParfor(max_clust);
% 
%     for n_cluster = 1:max_clust
%         send(q, n_cluster);
% 
%         idx_cluster = find(clusters == n_cluster);
% 
%         % Skip if only one mz inside cluster
%         if numel(unique(feature_table(idx_cluster,2))) < 2
%             continue
%         end
% 
%         % ---------------------------------------------------------
%         % Collapse cluster: keep max(volume) per SAMPLE
%         % ---------------------------------------------------------
%         % sample_vec_tmp = unique(feature_table(idx_cluster,5));
%         [~,~,u_ind] = unique(feature_table(idx_cluster,2));
%        %  count = accumarray(u_ind,feature_table(idx_cluster,6)>0);
%        %  [~,max_ind] = max(count);
%        %   feature_table(idx_cluster,1)
%        % not_max =  feature_table(idx_cluster,2) ~= u(max_ind);
%         [gc, grps] = groupcounts(feature_table(idx_cluster,2));
% 
%         [~, a] = max(gc);
%         if nnz(gc(a) == gc)>1
%              idx = find(gc(a) == gc);
%              sum_components_all_samples = accumarray(u_ind, feature_table(idx_cluster,6));
%              [~,a]  = max(sum_components_all_samples(idx));
%              a = idx(a);
%         end 
%         not_max = any(grps(~ismember(grps,grps(a)))' == feature_table(:,2),2);
% 
%         % ind_sample = grps(a) == feature_table(:,2);
%         feature_table(not_max,:) = 0;
%           if ~isempty(feature_groups_all)
%         %Feature group all struct put empty
%         fn = fieldnames(feature_groups_all);
%         empty_struct = cell2struct(cell(size(fn)), fn, 1);
%         feature_groups_all(not_max) = empty_struct;
%           end
% 
%     end
% 
%     % -------------------------------------------------------------
%     % Remove zero rows created during collapsing
%     % -------------------------------------------------------------
% 
%     ind_zeros = feature_table(:,1) == 0;
%     feature_table(:,[1,3:4]) = feature_table_old(:,[1,3:4]);
%     feature_table(:,2) = clusters;
% 
%     feature_table(ind_zeros,:) = [];
%       if ~isempty(feature_groups_all)
%     feature_groups_all(ind_zeros) = [];
%       end
%     % store old cluster values
%     clusters_all_samples_before_duplicate_removal = clusters_all_samples_before_duplicate_removal(~ind_zeros);
% 
%     %New cluster values
%     clusters_all_samples = feature_table(:,2);
% 
%     %Assign correct cluster number 
%     for n = 1:length(feature_groups_all)
%         feature_groups_all(n).sample_feature_ID = feature_table(n,2);
%     end 
% 
% end
function [feature_table, feature_groups_all, clusters_all_samples,clusters_all_samples_before_duplicate_removal] = cleanup_clusters_feature_table(feature_groups_all,clusters_all_samples,feature_table, Options)

  
    % Assign same Rt to all with the same cluster: 
    feature_table_old = feature_table;
    Rt_median(:,1) = accumarray(feature_table(:,2),feature_table(:,3),[],@median);
    Rt_median(:,2) = accumarray(feature_table(:,2),feature_table(:,3),[],@median);
    Rt_median = Rt_median(feature_table(:,2),:);
    
    mz_median = accumarray(feature_table(:,2),feature_table(:,1),[],@median);
    mz_median = mz_median(feature_table(:,2));
    
    % Sort table by mz, RT, sample
    feature_table(:,1) =  mz_median ;
    feature_table(:,3:4) = Rt_median(:,:);
    [feature_table, sort_idx] = sortrows(feature_table, [1,3,4], "ascend");
    clusters_all_samples_before_duplicate_removal = clusters_all_samples(sort_idx);
    feature_table_old = feature_table_old(sort_idx,:);
    if ~isempty(feature_groups_all)
    feature_groups_all = feature_groups_all(sort_idx);
    end
    
    % Rt_median = Rt_median(sort_idx,:);
    % mz_median = mz_median(sort_idx);

    N = size(feature_table,1);
    clusters = zeros(1,N);
    clusters(1) = 1;
    

    % -------------------------------------------------------------
    % Build clusters based on RT proximity
    % -------------------------------------------------------------

    m = 1;
    % wU = [(1 + Options.ppm_dev*1e-6) * feature_table(m,1), Options.Rtdev];
    wU = [0, Options.Clustering.distance_componentization];
    for n = 2:N
        if all(feature_table(n,1) - feature_table(m,1) == wU(1) & feature_table(n,3:4) - feature_table(m,3:4) <= wU(2:3)) % all(abs(feature_table(n,[1,3,4]) - feature_table(m,[1,3,4])) <= wU)
            continue
        end
        m = n;
        % max_dist_rt = Rt_median(m,:) + Options.Clustering.distance_componentization;
        % wU(1) = (1 + Options.ppm_dev*1e-6) * feature_table(m,1);
        clusters(1,n) = 1;
    end
    
    clusters = cumsum(clusters);

    % -------------------------------------------------------------
    % PROCESS clusters without deleting inside the loop
    % -------------------------------------------------------------
    max_clust = max(clusters);
  
    q = progressParfor(max_clust);

    for n_cluster = 1:max_clust
        send(q, n_cluster);

        idx_cluster = find(clusters == n_cluster);

        % Skip if only one mz inside cluster
        if numel(unique(feature_table(idx_cluster,2))) < 2
            continue
        end

        % ---------------------------------------------------------
        % Collapse cluster: keep max(volume) per SAMPLE
        % ---------------------------------------------------------
        % sample_vec_tmp = unique(feature_table(idx_cluster,5));
        [~,~,u_ind] = unique(feature_table(idx_cluster,2));
       %  count = accumarray(u_ind,feature_table(idx_cluster,6)>0);
       %  [~,max_ind] = max(count);
       %   feature_table(idx_cluster,1)
       % not_max =  feature_table(idx_cluster,2) ~= u(max_ind);
        [gc, grps] = groupcounts(feature_table(idx_cluster,2));
        
        [~, a] = max(gc);
        if nnz(gc(a) == gc)>1
             idx = find(gc(a) == gc);
             sum_components_all_samples = accumarray(u_ind, feature_table(idx_cluster,6));
             [~,a]  = max(sum_components_all_samples(idx));
             a = idx(a);
        end 
        not_max = any(grps(~ismember(grps,grps(a)))' == feature_table(:,2),2);
        
        % ind_sample = grps(a) == feature_table(:,2);
        feature_table(not_max,:) = 0;
          if ~isempty(feature_groups_all)
        %Feature group all struct put empty
        fn = fieldnames(feature_groups_all);
        empty_struct = cell2struct(cell(size(fn)), fn, 1);
        feature_groups_all(not_max) = empty_struct;
          end
       
    end

    % -------------------------------------------------------------
    % Remove zero rows created during collapsing
    % -------------------------------------------------------------
    
    ind_zeros = feature_table(:,1) == 0;
    feature_table(:,[1,3:4]) = feature_table_old(:,[1,3:4]);
    feature_table(:,2) = clusters;
    
    feature_table(ind_zeros,:) = [];
      if ~isempty(feature_groups_all)
    feature_groups_all(ind_zeros) = [];
      end
    % store old cluster values
    clusters_all_samples_before_duplicate_removal = clusters_all_samples_before_duplicate_removal(~ind_zeros);
    
    %New cluster values
    clusters_all_samples = feature_table(:,2);

    %Assign correct cluster number 
    for n = 1:length(feature_groups_all)
        feature_groups_all(n).sample_feature_ID = feature_table(n,2);
    end 

end