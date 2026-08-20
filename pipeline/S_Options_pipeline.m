%Define Options 
% cd C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection 
% c
% addpath C:\Users\mht541\Documents\2D_LC_IMS_R\2DIMS\Matlab\PeakCapacity
% addpath C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection\preprocessing;
%read IS data
IS_mz = readtable(meta_data_path,'Sheet',"IS");
Options.IS.name = IS_mz{:,1};
input_parameters = readtable(meta_data_path,'Sheet',"input_parameters");
mz_IS = IS_mz.x_M_H_;
mz_nonIS = IS_mz.NonDeuterated_Protonated;
get_ref_spec_IS

%%
% Initialize results as empty matrix
Options.min_peak_distance = input_parameters{1,"Options_value"};%3;
Options.int_thresh =  input_parameters{2,"Options_value"};%3000;
Options.ppm_dev =  input_parameters{3,"Options_value"};
Options.Adducts =  [input_parameters{:,"Adducts_value"}]';%[	18.034374132,	22.98976928  ,38.96370648];
Options.Adducts =  Options.Adducts(~isnan(Options.Adducts));%[	18.034374132,	22.98976928  ,38.96370648];
Options.doPreFilter = input_parameters{21,"Options_value"} == 1;
Options.doMFonRaw_data = input_parameters{22,"Options_value"} == 1;
Options.mz.C13      =  input_parameters{4,"Options_value"};
Options.mzHydrogen  =  input_parameters{5,"Options_value"};

% Clean up
Options.num_fragments = input_parameters{6,"Options_value"};

%Clustering
Options.Clustering.cutoff =  1-input_parameters{7,"Options_value"};
% Options.Clustering.max_peak_distance = 3;
Options.Clustering.distance_componentization = input_parameters{8:9,"Options_value"}';

Options.Rtdev = input_parameters{8:9,"Options_value"}';
Options.Rtdev_within_sample = input_parameters{10:11,"Options_value"}';
Options.doPlot = false;

Options.doSmoothing = input_parameters{12,"Options_value"} == 1;

% [peak_borders,sum_elution_profile,Z_feature,allPeaks,cosMat,cosMat_combined,clusters,mass_spectrum]  = deal(cell(length(MSroi_aug),1));
% mz_candidate = find(mzroi_aug >= input_parameters{13,"Options_value"} & mzroi_aug <= input_parameters{14,"Options_value"})';
%% Smoothing shape
F1_input = input_parameters{15:16,"Options_value"}';
F2_input = input_parameters{17:18,"Options_value"}';
F1 = gausswin(F1_input(1),F1_input(2));
F2 = gausswin(F2_input(1),F2_input(2));
Options.Filter = F1*F2';
Options.Filter = Options.Filter./max(Options.Filter,[],'all');%mexican_hat_2d_aniso(3,1, 9, 3);%

%% Curve resolution 
algorithm = {'als','mult'};
Options.curve_resolution.algorithm  = algorithm{input_parameters{1,"curve_resolutions_value"}};
Options.curve_resolution.max_factor = input_parameters{2,"curve_resolutions_value"};
Options.curve_resolution.exp_var_thresh = input_parameters{3,"curve_resolutions_value"};
%% Load Rt
% rt = load('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\Rts_2D.mat');



%% Folding of 2D plot
Options.PhaseShift = input_parameters{19,"Options_value"};
Options.modTime = input_parameters{20,"Options_value"};
%% save
Options.dir.results = [input_parameters{3,"Paths_dir"},test_name];

