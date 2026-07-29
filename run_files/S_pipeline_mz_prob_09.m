%% Load meta data
% clear 
cd C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection
ProjectPBA
addpath('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection')
addpath C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\tensor_toolbox


cd('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\SparseTensor\')
data_dir = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\ROI\Samples_aug\';
% load([data_dir,'MSroi_aug.mat'])
load([data_dir,'mzroi_aug.mat'])

matFile_list = dir('*.mat');
S_Options_peak_detection_mz_prob_09
n_samples = length(matFile_list);

% Build fileList with short names
ftemp = struct;
for k = 1:length(matFile_list)

    ftemp(k).name = matFile_list(k).name(1:9);
    ftemp(k).folder = Options.dir.results;
end

%% Do peak detection 
tic
clc
q = progressParfor(n_samples);
parfor k = 1:n_samples
    send(q, k);
    
    %Load the sample
    d = load([matFile_list(k).folder,'\',matFile_list(k).name], 'Z' );
    
    %Run pipeline for one sample
    run_peak_detection_pipeline(k, d.Z, mz_candidate, Options, mz_IS, mz_nonIS, full(mzroi_aug), mzTarget,ftemp, refSpec);

end
toc
%% Make strute to group samples
fprintf(1,'Build feature_groups_all\n');
feature_groups_all = build_points_all_samples(ftemp,Options,mzroi_aug);

%% %--- Run clustering using reference sample ---
fprintf(1,'Group components across samples\n');
clusters_all_samples = referenceClustering_all_samples(feature_groups_all,Options);
fprintf(1,'Save results\n');
save(['C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\',test_name,'.mat'], "clusters_all_samples", "feature_groups_all","-v7.3","-nocompression");
numClusters = max(clusters_all_samples);
disp(['Number of clusters found: ', num2str(numClusters)]);

%%
toc