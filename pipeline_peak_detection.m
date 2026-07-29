cd C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection
addpath C:\Users\mht541\Documents\2D_LC_IMS_R\2DIMS\Matlab\PeakCapacity
%read IS data

IS_mz = readtable("IS_table.xlsx");
mz_IS = IS_mz.x_M_H_;
mz_nonIS = IS_mz.NonDeuterated_Protonated;

Options.num_smooth_points = [3,5];

% Initialize results as empty matrix
Options.min_peak_distance = 3;
Options.int_thresh = 3000;
Options.ppm_dev = 10;
Options.Adducts = [	18.034374132,	22.98976928  ,38.96370648];
Options.mzHydrogen = 1.0078250319;

% Clean up
Options.num_fragments = 1;

Options.Clustering.cutoff =  0.3;
% Options.Clustering.max_peak_distance = 3;
Options.Clustering.distance_componentization = [4,15];

Options.Rtdev = [5,15];
Options.doPlot = true;


[peak_borders,sum_elution_profile,Z_feature,allPeaks,cosMat,cosMat_combined,clusters,mass_spectrum]  = deal(cell(length(MSroi_aug),1));
mzTarget =mz_IS(12);275.235;%268.191;

[~, a] = min(abs(mzroi_aug - mzTarget));
mz_candidate = a:a+10000%;1:length(mzroi_aug);
Options.dir.results = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\Single_samples';

0;
close all
clc
% Run peak detection
for k = 1;%sample_vec
    tic
%       fprintf(1,'Peak detection sample: %i/%i\n',k,length(sample_vec))

    [peak_borders{k},sum_elution_profile{k},Z_feature{k},allPeaks{k}] = run_peak_detection(Z,mz_candidate,Options);
    toc 
    tic 
    % Clustering
    if size(allPeaks{k},1)>1
        %Calculate cosine similarity it does it in a blocked fashion to
        %avoid OOM
        points = build_points_struct(allPeaks{k}, sum_elution_profile{k});

%         fprintf(1,'Clustering sample: %i/%i\n',k,length(sample_vec))
        clusters{k} = referenceClustering_one_sample(points,Options);
        
        %detect internal standards and give them a new cluster number as
        %well as flagging them in allPeaks{k}(ind,7)
        for nIS = 1:length(mz_IS)
            [clusters{k},allPeaks{k}] = split_clusters_wIS([mz_IS(nIS),mz_nonIS(nIS)], k, allPeaks, clusters, mzroi_aug, Options);
        end
        toc
        
%         fprintf(1,'Cleanup clusters sample: %i/%i\n',k,length(sample_vec))
        tic
        % Clean up after the loop
        [peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, ~, cosMat_combined] = ...
            cleanup_clusters(k, peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, [], []);
        mass_spectrum{k}=accumarray([allPeaks{k}(:,4),clusters{k}],allPeaks{k}(:,5),[length(mzroi_aug),max(clusters{k})],[],[],true);
    else
        clusters{k} = 1;
    end

    %Filter out feature groups with <= Options.num_fragments
    [peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, ~, cosMat_combined] = ...
        filterClustersAndVars(peak_borders, sum_elution_profile, Z_feature, allPeaks, clusters, [], cosMat_combined, Options);

    % Visualize
    if Options.doPlot
        figure
        [~,fG]=min(abs(mzroi_aug(allPeaks{k}(:,4))-mzTarget));
        vizualize_feature_group(Z_feature{k},peak_borders{k},mzroi_aug,allPeaks{k},clusters{k},fG)
    end
    
%     save_sample_results(fileList(k).name(1:end-10),Options,peak_borders{k}, sum_elution_profile{k}, Z_feature{k}, allPeaks{k}, clusters{k}, cosMat_combined,mzroi_aug)
   clc
   toc 
end
