%Define Options 
cd C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection
addpath C:\Users\mht541\Documents\2D_LC_IMS_R\2DIMS\Matlab\PeakCapacity
addpath C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection\preprocessing;
%read IS data

IS_mz = readtable("IS_table.xlsx");
mz_IS = IS_mz.x_M_H_;
mz_nonIS = IS_mz.NonDeuterated_Protonated;
% get_ref_spec_IS


% Initialize results as empty matrix
Options.min_peak_distance = 3;
Options.int_thresh = 3000;
Options.ppm_dev = 10;
Options.Adducts = [	18.034374132,	22.98976928  ,38.96370648];
Options.mzHydrogen = 1.0078250319;

% Clean up
Options.num_fragments = 0;

Options.Clustering.cutoff =  0.3;
% Options.Clustering.max_peak_distance = 3;
Options.Clustering.distance_componentization = [4,15];

Options.Rtdev = [5,15];
Options.doPlot = false;

Options.doSmoothing = true;

% [peak_borders,sum_elution_profile,Z_feature,allPeaks,cosMat,cosMat_combined,clusters,mass_spectrum]  = deal(cell(length(MSroi_aug),1));
mzTarget = mz_IS(6);275.235;%268.191;
[~, a] = min(abs(mzroi_aug - mzTarget));

mz_candidate = [1:length(mzroi_aug)];
%%
F1 = gausswin(7,6);
F2 = gausswin(19,4);
Options.Filter = F1*F2';
Options.Filter = Options.Filter./sum(Options.Filter,'all');%mexican_hat_2d_aniso(3,1, 9, 3);%
% mesh(Options.Filter)
%%
test_name = 'No_fragmentFilter_2DSmoothing_25927';
Options.dir.results = ['C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\',test_name];
