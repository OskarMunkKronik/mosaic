clf

% feature_list = clusters_05(6).names(1);2933;
load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\curve_resolution_final_meta_data_Ver2.mat")
peak_table = readtable("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\peak_table_final.xlsx","Sheet","peak_volume",'Range','A3:C18005');
load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\curve_resolution_final_251110_Ver2.mat","clusters_all_samples_before_duplicate_removal","feature_table")
load_dir = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\component_cluster\';
n_samples = 78;
cd C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection
ProjectPBA
metadata = readtable("H:\PhD\Articles\8_LCxLC_Vandalf_BigDataSet\Part_2_Data\Kristoffer_data\Kristoffer_ID_catchment.xlsx",'ReadVariableNames',true)
%%
clc
clf
figure(1)
col = lines(n_samples);
mzTarget = 271.0318




% mz_IS(1)
feature_list = find((abs(mzTarget-peak_table{:,3})./peak_table{:,3}).*1e6 < Options.ppm_dev);

plot_elution_profile_and_mass_spectrum( ...
    feature_list, ...
    feature_table, ...
    clusters_all_samples_before_duplicate_removal, ...
    mzroi_aug, ...
    load_dir, ...
    Options, ...
    n_samples, ...
    col, ...
    1:78);
