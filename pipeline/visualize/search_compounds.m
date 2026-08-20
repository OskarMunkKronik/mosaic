%% ================== SETUP ==================
dataDir = Options.Paths.save2mat;
load(fullfile(dataDir,'curve_resolution_final_meta_data_Ver2.mat'),'Options')
load(fullfile(Options.Paths.save2mat,'feature_table_before_duplicate_remove.mat'))
load(fullfile(Options.Paths.save2mat,'curve_resolution_final.mat'),'clusters_all_samples_before_duplicate_removal')
load(fullfile(Options.Paths.save2mat,'ROI\Samples_aug\MSroi_aug.mat'),"Rt_aug")
load(fullfile(Options.Paths.save2mat,'ROI\Samples_aug\mzroi_aug.mat'))
load(fullfile(Options.Paths.save2mat,"ROI\fileList_Samples.mat"))
% load(fullfile(Options.Paths.save2mat,'curve_resolution_final_meta_data_Ver2.mat'))
[clusters_all_samples, feature_groups_all] = load_MF_results_parts( fullfile(Options.Paths.save2mat,'NMF_final'),'NMF_final_part_')
% load(fullfile(Options.Paths.save2mat,'mass_spectra.mat'))
load_dir = fullfile(Options.Paths.save2mat,'component_cluster');

n_samples = max(feature_table(:,5));
num_clusters_before_curve_resolution = length(mz_ind);
% Add utils path (safer than cd)
% cd('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection')
% ProjectPBA

%% ================== MASS SPECTRA PREP ==================
mass_spectra_cell = {feature_groups_all.massSpec};

% Remove zeros and map mz
valid_idx = cellfun(@(s) s > 0, mass_spectra_cell(:), 'UniformOutput', false);

mz_in_mass_spectra = cellfun(@(mz,idx) mz(idx), ...
    repmat({mzroi_aug}, size(valid_idx)), valid_idx, 'UniformOutput', false);

mass_spectra_cell = cellfun(@(spec,idx) spec(idx), ...
    mass_spectra_cell, valid_idx, 'UniformOutput', false);

len_spectra = cellfun(@numel, mz_in_mass_spectra);

%% ================== TARGETS ==================
ppm_dev = 50; % Set the ppm deviation to a user defined value corresponding to the mass accuracy of your mass spectrometer
sample_names = {fileList.name}';
rt_offset = round(median(cellfun(@(rt) rt(1),Rt_aug))/60,3);

%% ================== Search ==================

compound_name = {'non'}; % Fill in for plotting 
% Plot components with mzTarget in their Mass spectrum
mzTarget =447.9899





; % Target m/z to search mass spectra for 
figure(1)
%figure
search_plot_volume_elu_ms( ...
    num_clusters_before_curve_resolution, mz_in_mass_spectra, mzTarget, ppm_dev, ...
    feature_table, rt_offset, Options, ...
    clusters_all_samples_before_duplicate_removal, ...
    mzroi_aug, load_dir, n_samples, ...
    sample_names, compound_name);

%% Plot specific components
figure(2)
42
componentTarget = 1624% The specific component desired to be plotted
volume = plot_volume_elu_ms_component( ...
    num_clusters_before_curve_resolution,componentTarget, ...                         % <-- directly provided
    feature_table, rt_offset, Options, ...
    clusters_all_samples_before_duplicate_removal, ...
    mzroi_aug, load_dir, n_samples, ...
    sample_names, []);
