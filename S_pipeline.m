clear
addpath('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection')
addpath C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\tensor_toolbox

% do peak detection
cd('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\SparseTensor\')
data_dir = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\ROI\Samples_aug\';
% load([data_dir,'MSroi_aug.mat'])
load([data_dir,'mzroi_aug.mat'])

matFile_list = dir('*.mat');
S_Options_peak_detection
n_samples = length(matFile_list);

%% Do peak detection 
tic
for k = 1:n_samples
    load([matFile_list(k).folder,'\',matFile_list(k).name])
    S_Options_peak_detection
    [~,a] = min(abs(mzroi_aug-mzTarget));
       run_peak_detection_pipeline(k, Z, mz_candidate, Options, mz_IS, mz_nonIS, mzroi_aug, mzTarget, fileList  );

end
%% Make strute to group samples
ftemp = struct;
for k = 1:length(matFile_list)
    ftemp(k).name = matFile_list(k).name(1:9);
    ftemp(k).folder = Options.dir.results;
end
feature_groups_all = build_points_all_samples(ftemp,Options,mzroi_aug);

%% %--- Run clustering using reference sample ---
clusters_all_samples = referenceClustering_all_samples(feature_groups_all,Options);
save('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\Feature_groups_all_samples_No_Fragment_Filter_No_smoothing.mat', 'clusters_all_samples', 'feature_groups_all');
numClusters = max(clusters_all_samples);
disp(['Number of clusters found: ', num2str(numClusters)]);

%%
toc