%% ================== SETUP ==================
dataDir = 'D:\D4RUNOFF\Mosaic_300';%Options.Paths.save2mat;
load(fullfile(dataDir,'curve_resolution_final_meta_data_Ver2.mat'),'Options')
load(fullfile(Options.Paths.save2mat,'feature_table_before_duplicate_remove.mat'))
load(fullfile(Options.Paths.save2mat,'curve_resolution_final.mat'),'clusters_all_samples_before_duplicate_removal')
load(fullfile(Options.Paths.save2mat,'ROI\Samples_aug\MSroi_aug.mat'),"Rt_aug")
load(fullfile(Options.Paths.save2mat,'ROI\Samples_aug\mzroi_aug.mat'))
load(fullfile(Options.Paths.save2mat,"ROI\fileList_Samples.mat"))
load(fullfile(Options.Paths.save2mat,'curve_resolution_final_meta_data_Ver2.mat'))
load(fullfile(Options.Paths.save2mat,'mass_spectra.mat'))
load_dir = fullfile(Options.Paths.save2mat,'component_cluster');

n_samples = max(feature_table(:,5));
num_clusters_before_curve_resolution = length(mz_ind);
% Add utils path (safer than cd)
% cd('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection')
% ProjectPBA

%% ================== MASS SPECTRA PREP ==================
mass_spectra_cell = {feature_groups_all.massSpec};

% Remove zeros and map mz
valid_idx = cellfun(@(s) s > 0, mass_spectra_cell, 'UniformOutput', false);

mz_in_mass_spectra = cellfun(@(mz,idx) mz(idx), ...
    repmat({mzroi_aug}, size(valid_idx)), valid_idx, 'UniformOutput', false);

mass_spectra_cell = cellfun(@(spec,idx) spec(idx), ...
    mass_spectra_cell, valid_idx, 'UniformOutput', false);

len_spectra = cellfun(@numel, mz_in_mass_spectra);

%% ================== TARGETS ==================

ppm_dev = 200;
sample_names = {fileList.name}';
% sample_names = cellfun(@(c) c(32:end-11), sample_names, 'UniformOutput', false);
rt_offset = round(median(cellfun(@(rt) rt(1),Rt_aug))/60,3);
sample_names = {'Standards BioOils',...
'P#6',...
'Feed 1',...
'P#11',...
'QC BTG all 20 ppm',...
'KU1'}
% sample_names = contains(sample_names,"20250916_Test_Aragorn~20250916-")
%% ================== Search ==================

compound_name = {'non'}
% Plot components with mzTarget in their Mass spectrum
mzTarget =   462.0071
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
componentTarget = 1023;

volume = plot_volume_elu_ms_component( ...
    num_clusters_before_curve_resolution,componentTarget, ...                         % <-- directly provided
    feature_table, rt_offset, Options, ...
    clusters_all_samples_before_duplicate_removal, ...
    mzroi_aug, load_dir, n_samples, ...
    sample_names, []);
% set(gcf, 'Color', 'w')   % Set figure background to white
% set(gca, 'Color', 'w')   % Set axes background to white