% All pipeline
meta_data_path = "D:\Nadine\Mosaic_300\meta_data_300_counts.xlsx"
% meta_data_path = "F:\Nadine\meta_data"%thresh_morethan1mz.xlsx"
%% Get all paths
fprintf(1,'Setup Mosaic\n');
cd C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection
ProjectPBA

%% run options file
test_name = ['GCxGC_300'];
S_Options_pipeline
S_OptionsStruct_ROI

% addpaths 
addpath(input_parameters{6,"Paths_dir"}{:})
%% ROI
fprintf(1,'ROI detection\n');
S_Import_template_ROI

%% Fold to tensor
fprintf(1,'Convert to sparse tensor\n');
S_Convert2spTensor

%% Run MF 
fprintf(1,'Run Mass filtering\n');
S_pipeline_feature_detection_MF

%% Run curve-resolution
fprintf(1,'Run curve-resolution on Mass filtered data\n');
S_pipeline_run_curve_resolution_all_samples

%% Export files 
fprintf(1,'Export files...\n');
S_pipeline_export_feature_table_to_excell
S_save_mass_spectra