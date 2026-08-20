% All pipeline
meta_data_path = "N:\SCIENCE-PLEN-ECP_AnalytChem\People\Nadine\publication\meta_data.xlsx";
% meta_data_path = "F:\Nadine\meta_data"%thresh_morethan1mz.xlsx"
%% Get all paths
fprintf(1,'Setup Mosaic\n');
cd D:\people\mht541\projects\mosaic\mosaic_published
ProjectPBA

%% run options file
test_name = [''];
S_Options_pipeline
S_OptionsStruct_ROI

% addpaths 
addpath(input_parameters{6,"Paths_dir"}{:})
%% ROI
fprintf(1,'ROI detection\n');
cd([Options.Paths.CDF])
fileList                                          = dir(['*',Options.ROI.MS1_Suffix]);
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