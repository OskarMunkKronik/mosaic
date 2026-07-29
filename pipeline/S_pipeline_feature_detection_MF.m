%% Load meta data
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
tic
% clc
q = progressParfor(n_samples);
parfor k = 1:n_samples
    send(q, k);
    
    %Load the sample
    d = load([matFile_list(k).folder,'\',matFile_list(k).name], 'Z' );
    
    %Run pipeline for one sample
    run_peak_detection_pipeline(k, d.Z, mz_candidate, Options, mz_IS, mz_nonIS, full(mzroi_aug), [],ftemp, refSpec);

end
toc
%% Make strute to group samples
fprintf(1,'Build feature_groups_all\n');
feature_groups_all = build_points_all_samples(ftemp,Options,mzroi_aug);


%% %--- Run clustering using reference sample ---
tic
clc
fprintf(1,'Group components across samples\n');
clusters_all_samples = referenceClustering_all_samples(feature_groups_all,Options);
process_time.MF_process_time = toc;
%%
tic
fprintf(1,'Save results\n');
save(fullfile(Options.Paths.save2mat,['MF_results.mat']), "clusters_all_samples", "feature_groups_all","-v7.3","-nocompression");
numClusters = max(clusters_all_samples);
disp(['Number of clusters found: ', num2str(numClusters)]);
process_time.MF_process_time = toc;
%%