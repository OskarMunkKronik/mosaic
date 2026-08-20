%% Load meta data
Options.processing_time.MF.process.all = tic;
cd(fullfile(Options.Paths.save2mat,'SparseTensor'))
matFile_list = dir('*.mat');
% S_Options_peak_detection_mz_prob
n_samples = length(matFile_list);

% Build fileList with short names
ftemp = struct;
for k = 1:length(matFile_list)

    ftemp(k).name = matFile_list(k).name;
    ftemp(k).folder = Options.dir.results;
end

%% Do peak detection 
% tic
Options.processing_time.MF.process.peak_detection = tic;
% clc
q = progressParfor(n_samples);
parfor k = 1:n_samples
    send(q, k);
    
    %Load the sample
    d = load([matFile_list(k).folder,'\',matFile_list(k).name], 'Z' );
    
    %Run pipeline for one sample
    run_peak_detection_pipeline(k, d.Z, mz_candidate, Options, mz_IS, mz_nonIS, full(mzroi_aug), [],ftemp, refSpec);

end
Options.processing_time.MF.process.peak_detection = toc(Options.processing_time.MF.process.peak_detection);
%% Make strute to group samples
Options.processing_time.MF.process.componentization = tic;
fprintf(1,'Build feature_groups_all\n');
feature_groups_all = build_points_all_samples(ftemp,Options,mzroi_aug);
%% %--- Run clustering using reference sample ---
tic
clc
fprintf(1,'Group components across samples\n');
clusters_all_samples = referenceClustering_all_samples(feature_groups_all,Options);
feature_table = full([feature_groups_all.basePeak; clusters_all_samples'; feature_groups_all.Rt;feature_groups_all.sample;feature_groups_all.volume_summed]');

Options.processing_time.MF.process.componentization = toc(Options.processing_time.MF.process.componentization);
%% Save
Options.processing_time.MF.save = tic;
fprintf(1,'Save results\n');
% save(fullfile(Options.Paths.save2mat,['MF_results.mat']), "clusters_all_samples", "feature_groups_all","-v7.3");
save_MF_results_parts(clusters_all_samples, feature_groups_all, Options.Paths.save2mat,'MF_results_part', 10)
numClusters = max(clusters_all_samples);
disp(['Number of clusters found: ', num2str(numClusters)]);
Options.processing_time.MF.save =  toc(Options.processing_time.MF.save);
%%
Options.processing_time.MF.process.all = toc(Options.processing_time.MF.process.all);