
% %%
% close all
% tic
% vol_curve = zeros(n_samples,length(mz_IS));
% mz_vec = mz_IS;%full(unique([feature_groups_all.basePeak]))';%;
% % [~,refs] = max([feature_groups_all.volume_summed]);
% % i_vec = find(refs ~= 1:length(mz_vec));
% i_vec = 8;1:length(mz_vec);
%
% cluster_curve_resolution = zeros(length([feature_groups_all.volume_summed]),1);
% counter = 1;
% Results = struct();
% Options.curve_resolution.max_factor =  4;
% Options.curve_resolution.exp_var_thresh = 0.95;
% d = full([feature_groups_all.basePeak; clusters_all_samples'; feature_groups_all.Rt;feature_groups_all.sample;feature_groups_all.volume_summed;feature_groups_all.volume_basePeak]');
% for i = i_vec
%     % Find mz's, peak b
%     mzTarget = mz_vec(i);
%       if i == 1
%     mzTarget = mz_vec(i)+Options.Adducts(2)-Options.mzHydrogen;%275.2350%268.191
%     else
%     mzTarget = mz_vec(i);%275.2350%268.191
%     end
%
%     [mz_ind,peak_borders,gc,grps,Rt_cluster] = get_mz_ind_for_curve_resolution(d,mzTarget,clusters_all_samples,feature_groups_all,Options);
%
%     %% Select cluster closest in Rt to the IS
%     n_cluster = find(cellfun(@(b) ...
%         all(abs([ ...
%         Options.coords.X(1, round(mean(b(1,:),2))), ...
%         Options.coords.Y(round(mean(b(2,:),2)), 1) ...
%         ] - rt_exp(i,:)) < [0.9 1] ), ...
%         Rt_cluster));
%     if isempty(n_cluster)
%         n_cluster = 1;
%     end
%
%     % if n_cluster
%     %%
%     matFile_list_tmp = dir([Options.dir.results,'\*.mat']);
%     q = progressParfor(n_samples);
%     Z = zeros([diff(peak_borders{n_cluster},[],2)'+1,length(mz_ind{n_cluster}),n_samples]); %1D x 2D x mz x samples
%     for n_cluster = 1:length(mz_ind)
%     parfor k = 1:n_samples
%         send(q, k);
%
%         %Load the sample
%         data = load([matFile_list(k).folder,'\',matFile_list(k).name], 'Z' );
%         Z_filtered = double(data.Z(peak_borders{n_cluster}(1,1):peak_borders{n_cluster}(1,2), ...
%             peak_borders{n_cluster}(2,1):peak_borders{n_cluster}(2,2), ...
%             mz_ind{n_cluster}));
%         for i = 1:length(mz_ind{n_cluster})
%
%             Z_filtered(:,:,i) = conv2(Z_filtered(:,:,i), Options.Filter, 'same');
%         end
%         Z(:,:,:,k) = Z_filtered;
%     end
%     toc
%     %% Do curve-resolution
%
%     [W, H, cg, exp_var, X_unfolded] = perform_nnmf_unfold(Z, Options);
%     X_resolved = W*H;
%     figure
%     col = 'rbgmk';lines(max(cg));
%     [X,Y] = meshgrid(peak_borders{n_cluster}(1,1):peak_borders{n_cluster}(1,2), ...
%         peak_borders{n_cluster}(2,1):peak_borders{n_cluster}(2,2));
%     clear p
%     Z_refolded = cell(max(cg),1);
%     for n = 1:max(cg)
%
%         subplot(1,2,1)
%         [d1, d2, d3, d4] = size(Z);
%
%         Z_refolded{n} = zeros(d1, d2, d3, d4);
%
%         for nn = 1:d3
%             % Take column n and reshape it back to the permuted form
%             X_tmp = reshape(X_resolved(:,nn), [d2, d1, d4]);  % note: d2 first because of permute [2,1,3]
%             % Undo the permute to restore original order
%             Z_refolded{n}(:,:,nn,:) = permute(X_tmp, [2,1,3]);
%         end
%
%
%         h = 1;
%         % clear p
%         % for n_mz = 1:sum(cg == n)
%         for k = 1:n_samples
%             % Z_tmp_transpose = Z_tmp(:,:,n_mz , k);
%
%             p_tmp =  plot3(X,Y,squeeze(sum(Z_refolded{n}(:,:,cg == n,k),3)),'Color',col(n));
%
%             hold on
%         end
%         p(h) = p_tmp(1);
%
%         h = h +1;
%         % end
%         hold on
%         % legend(p,{num2str(mzroi_aug(mz_ind(cg == n)))})
%
%     end
%
%     ind = false(max(cg),1);
%     for n = 1:max(cg)
%
%         subplot(max(cg),2,n*2)
%         stem(mzroi_aug(mz_ind{n_cluster}(cg == n)),sum(X_unfolded(:,(cg == n))),'Color',col(n))
%         xlim([min(mzroi_aug(mz_ind{n_cluster})-10) max(mzroi_aug(mz_ind{n_cluster})+10)])
%         ppmError = abs(mzroi_aug(mz_ind{n_cluster}(cg == n)) - mzTarget) ./ mzTarget * 1e6;
%         ind(n) = nnz(ppmError < Options.ppm_dev);
%
%
%     end
%
%     vol_curve(:,i) = squeeze(sum(Z_refolded{n}(:,:,cg == find(ind,1,'first'),:),[1:3]));
%
%     sgtitle([IS_mz{i,1}{:},', mz: ', num2str(mzTarget)])
%
%     %% Build
%     % Fill struct entry
%         max_cluster = max(cluster_curve_resolution);
%     for n = 1:length(Z_refolded)
%         vol_tmp = squeeze(sum(Z_refolded{n}(:,:,cg == n,:),1:3));
%         for k = 1:n_samples
%             a = find(cg == n);
%             mass_spec = squeeze(sum(Z_refolded{n}(:,:,(cg == n),k),[1,2]));
%             [~,basepeak_ind] = max(mass_spec);
%
%             BPC = Z_refolded{n}(:,:,basepeak_ind,k);
%             Results(counter).feature_elu_2D    = BPC;
%             if nnz(BPC) > 0
%
%
%                 [r,c] = find(max(BPC,[],'all') == BPC);
%                 Results(counter).Rt                = [X(1,r);Y(c,1)];
%                 Results(counter).massSpec          = sparse(zeros(length(mzroi_aug),1));
%                 Results(counter).massSpec(mz_ind{n_cluster}(cg == n)) = sparse(mass_spec);
%                 Results(counter).volume_summed     = vol_tmp(k);
%                 Results(counter).volume_basePeak   = sum(BPC,'all');
%                 Results(counter).basePeak          = mzroi_aug(mz_ind{n_cluster}(a(basepeak_ind)));
%                 Results(counter).sample_name            = matFile_list_tmp(k).name;
%                 Results(counter).sample            = k;
%                 Results(counter).sample_feature_ID = counter;
%                 Results(counter).sample_cluster_ID = unique(cg(cg == n));
%                 cluster_curve_resolution(counter)  = Results(counter).sample_cluster_ID + max_cluster;
%                 counter = counter + 1;
%             end
%         end
%
%     end
%     end
% end
% if isempty(Results(counter).Rt)
%     Results(counter) = [];
% end
% cluster_curve_resolution = cluster_curve_resolution(1:counter-1);
%
%


