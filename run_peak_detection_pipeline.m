function [peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, mass_spectrum, cosMat_combined] = ...
    run_peak_detection_pipeline(k, Z, mz_candidate, Options, mz_IS, mz_nonIS, mzroi_aug, mzTarget, fileList , refSpec )


cosMat_combined =  [];
[peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters] = deal(cell(length(k),1));
% Run peak detection
%     fprintf(1,'Peak detection sample: %i\n', k);
[peak_borders{k}, sum_elution_profile{k}, Z_feature{k}, allPeaks{k}] = ...
    run_peak_detection(Z, mz_candidate, Options);

% --- Clustering ---
if size(allPeaks{k},1) > 1
    points = build_points_struct(allPeaks{k}, sum_elution_profile{k});
    %         fprintf(1,'Componentizing sample: %i\n', k);
    clusters{k} = referenceClustering_one_sample(points,mzroi_aug(allPeaks{k}(:,4)),Z_feature{k}, Options);


    % Detect internal standards
   if ~isempty(mz_IS)
    for nIS = 1:length(mz_IS)
        % [~,a] = min(abs(mzroi_aug - mz_IS(nIS)));
        % if a < allPeaks{k}(end,4)
            [clusters{k}, allPeaks{k}] = split_clusters_wIS( ...
                [mz_IS(nIS), mz_nonIS(nIS)], k, allPeaks, clusters, mzroi_aug, Options,refSpec{nIS});
        % end
    end
   end 

   
    % Cleanup
    %         fprintf(1,'Clean up components : %i\n', k);
   [peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, ~, cosMat_combined] = ...
        cleanup_clusters(k, peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, [], []);
    
    % Perform nnmf on peaks where the peak apex is not identical
    points = build_points_struct(allPeaks{k}, sum_elution_profile{k});
    Z_tmp = get_elution_profiles_cluster(Z_feature, peak_borders, clusters, points, k);
    flag = find(~cellfun(@isempty, Z_tmp));

    for n = 1:length(flag)
        [~, H, cg, exp_var, ~] = perform_nnmf_unfold(Z_tmp{flag(n)}, Options);
        cluster_members = find(clusters{k} == flag(n));
        clusters{k}(cluster_members(cg>1) ) = max(clusters{k}) + 1;
    end
     % Cleanup
    [peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, ~, cosMat_combined] = ...
        cleanup_clusters(k, peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, [], []);

    % Build mass spectrum
    mass_spectrum = accumarray([allPeaks{k}(:,4), clusters{k}], ...
        allPeaks{k}(:,5), [length(mzroi_aug), max(clusters{k})], [], [], true);
else
    clusters{k} = 1;
    mass_spectrum = [];
end

% --- Filter clusters ---
[peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, ~, cosMat_combined] = ...
    filterClustersAndVars(peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, [], cosMat_combined, Options);

% --- Visualization ---
% if  ~isempty(allPeaks{k})
%     if ~isempty(mz_IS)
%     close all
%     IS_name = Options.IS.name;
%     for nIS = 1:14
%         if exist(fullfile(Options.Paths.save2mat,'IS',IS_name{nIS}),'dir') ~= 7
%             mkdir(fullfile(Options.Paths.save2mat,'IS',IS_name{nIS}))
%         end 
%         mzTarget = mz_IS(nIS);
%         figure
%         [val, fG] = min(abs(mzroi_aug(allPeaks{k}(:,4)) - mzTarget));
%         if val/mzTarget *1e6 < Options.ppm_dev
%         vizualize_feature_group(Z_feature{k}, peak_borders{k}, mzroi_aug, allPeaks{k}, clusters{k}, fG);
%         sgtitle(mzTarget)
%         fig = gcf;
%         saveas(fig, fullfile(Options.Paths.save2mat,'IS',IS_name{nIS},[fileList(k).name,'.jpeg']));
%         end 
%     end
% 
% end
% end 
% --- Save results ---
save_sample_results(fileList(k).name, Options, ...
    peak_borders{k}, sum_elution_profile{k}, Z_feature{k}, allPeaks{k}, clusters{k}, cosMat_combined, mzroi_aug, mass_spectrum);

end
