% All pipeline
clear Options
meta_data_path = "C:\Users\mht541\Desktop\test_mosaic\pipeline\meta_data.xlsx";
% meta_data_path = "F:\Nadine\meta_data"%thresh_morethan1mz.xlsx"
%% Get all paths
fprintf(1,'Setup Mosaic\n');
cd C:\Users\mht541\Desktop\test_mosaic
ProjectPBA

%% run options file
test_name = [''];


for n_samples_speed_test = [1,5,15,30,50,78]
    cd C:\Users\mht541\Desktop\public_mosaic
    S_Options_pipeline
S_OptionsStruct_ROI

% addpaths 
addpath(input_parameters{6,"Paths_dir"}{:})
%% ROI
fprintf(1,'ROI detection\n');
    paths = {
    'C:\Users\mht541\Desktop\public_mosaic\curve_resolution_results'
    'C:\Users\mht541\Desktop\public_mosaic\MF'
    'C:\Users\mht541\Desktop\public_mosaic\ROI'
    'C:\Users\mht541\Desktop\public_mosaic\SparseTensor'
    'C:\Users\mht541\Desktop\public_mosaic\TIC'
    'C:\Users\mht541\Desktop\public_mosaic\curve_resolution_final.mat'
    'C:\Users\mht541\Desktop\public_mosaic\curve_resolution_final_meta_data_Ver2.mat'
    'C:\Users\mht541\Desktop\public_mosaic\feature_groups_all_before_duplicate_remove.mat'
    'C:\Users\mht541\Desktop\public_mosaic\feature_table_before_duplicate_remove.mat'
    'C:\Users\mht541\Desktop\public_mosaic\mass_spectra.mat'
    'C:\Users\mht541\Desktop\public_mosaic\mass_spectra.msp'
    'C:\Users\mht541\Desktop\public_mosaic\MF_results.mat'
    'C:\Users\mht541\Desktop\public_mosaic\BPC'
    'C:\Users\mht541\Desktop\public_mosaic\clusters_for_curve_resolution'
    'C:\Users\mht541\Desktop\public_mosaic\component_cluster'
};

for i = 1:numel(paths)
    p = paths{i};

    if isfolder(p)
        rmdir(p, 's');   % Delete folder and all its contents
        fprintf('Deleted folder: %s\n', p);
    elseif isfile(p)
        delete(p);       % Delete file
        fprintf('Deleted file: %s\n', p);
    else
        fprintf('Not found: %s\n', p);
    end
end
    Options.processing_time.mosaic.process.all = tic;
cd([Options.Paths.CDF])
fileList                                          = dir(['*',Options.ROI.MS1_Suffix]);
fileList = fileList(1:n_samples_speed_test);%randi([1,length(fileList)],n,1));
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
Options.processing_time.mosaic.process.all = toc(Options.processing_time.mosaic.process.all);
save(fullfile(Options.Paths.save2mat,[ 'Options_nSamples_',num2str(n_samples_speed_test),'.mat']),"fileList","Options")
plot_processing_times
end 