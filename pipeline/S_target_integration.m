% All pipeline
meta_data_path ="N:\SCIENCE-PLEN-ECP_AnalytChem\Projects\MST\mosaic\Batch_2_influent\260923_meta_data_nominal_mz.xlsx";
% meta_data_path = "F:\Nadine\meta_data"%thresh_morethan1mz.xlsx"
%% Load input data
input_data = readtable("N:\SCIENCE-PLEN-ECP_AnalytChem\Projects\MST\mosaic\target_integration_template.xlsx")
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

cd C:\Users\mht541\Documents\MATLAB
delete(gcp('nocreate'))
parpool(6)
S_Import_template_ROI
cd C:\Users\mht541\Documents\MATLAB
delete(gcp('nocreate'))
parpool(6)
%% Fold to tensor
fprintf(1,'Convert to sparse tensor\n');
S_Convert2spTensor

%% Integration
rt_expected = input_data{:,["x1stDimensionTime_min_","x2ndDimensionTime_s_"]};
rt_tol = input_data{:,["x1stDimensionR_T_Deviation_s_","x2ndDimensionR_T_Deviation_s_"]}./[60,1];
mzTarget = input_data{:,"mzTarget"};
tol = input_data{:,"mzTol"}; 
compound_name = input_data{:,"Name"};

[volume_manual, apexTables] = extractTargetVolumes(Options, fileList, ...
        mzTarget, tol, rt_expected, rt_tol, compound_name,1:length(fileList));

%% Reprocessing volumes if rt borders need to change depending on co-elution. 
volume_manual = reprocessVolumesFromExcel(Options, matFile_list, tol, 5);
