cd C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection
ProjectPBA

cd C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection\run_files\

data_dir = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\ROI\Samples_aug\';
load([data_dir,'mzroi_aug.mat'])
S_Options_peak_detection_mz_prob
load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\curve_resolution_final_251110_Ver2.mat","feature_table")
load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\curve_resolution_final_meta_data_Ver2.mat")

%% --- Load Metadata ---
SampleList = readtable("H:\PhD\Experiments\LargeSampleSet_Vandalf\SampleList\SampleList_Combined.xlsx");
fileList_tmp = matFile_list;

% Keep only matching samples
SampleList = SampleList(contains(SampleList.SampleName_1, extractBefore({fileList_tmp.name},10)),:);

% Remove specific sample if needed
SampleList = SampleList(~contains(SampleList.SampleName_1,{'270724_38.d'}),:);

% Define groups
vec_QC     = find(ismember(SampleList.WWTP,{'QC1','QC2','QC3'}));
vec_ValQC  = find(ismember(SampleList.WWTP,{'ValQC1','ValQC2'}));
vec_blank  = find(ismember(SampleList.WWTP,{'Blank'}));
vec_Source = find(ismember(SampleList.WWTP,{'SA','SB'}));

% Exclude groups from fileList
exclude_idx = [vec_blank; vec_QC; vec_Source; vec_ValQC];
fileList_tmp(exclude_idx) = [];

% Indices for PCA samples
sample_vec_pca = setdiff(1:size(SampleList,1), exclude_idx);

WWTP_type = unique(SampleList.WWTP);
sample_id = SampleList.WWTP;
%% Calculate for table 
T = accumarray([feature_table(:,2),feature_table(:,5)],feature_table(:,6));

Rt(:,1) = accumarray(feature_table(:,2),Options.coords.X(1,feature_table(:,3)),[],@median);
Rt(:,2) = accumarray(feature_table(:,2),Options.coords.Y(feature_table(:,4),1),[],@median);
basepeak = accumarray(feature_table(:,2),feature_table(:,1),[],@median);
Component_group = accumarray(feature_table(:,2),feature_table(:,2),[],@median);
detection_rate_all_samples =round(sum(T>0,2)*100./size(T,2));
idx_QC = contains(sample_id,'QC') & ~contains(sample_id,'ValQC') ;
QC_std = std(T(:,idx_QC),[],2);
QC_mean = mean(T(:,idx_QC),2);
rsd_qc = QC_std./QC_mean .* 100;
blank_mean = mean(T(:,contains(sample_id,'Blank') & ~contains(sample_id,'InstrumentBlank')),2);
max_int = max(T(:,~contains(sample_id,{'QC','Blank'})) ,[],2);
min_int = min(T(:,~contains(sample_id,{'QC','Blank'})) ,[],2);
% Combine into one matrix

%% Export to table 
clc
% numeric part
M_numeric = [Component_group, detection_rate_all_samples, basepeak, Rt,QC_std,QC_mean,rsd_qc,blank_mean,min_int,max_int, T];

% header row as strings
header = [{'#', '%', 'm/z', '¹tᵣ (min)', '²tᵣ (sec)',['summed peak volume'],['summed peak volume'],['%'],'summed peak volume','summed peak volume','summed peak volume'}, sample_id{:}];

% convert numeric matrix to a cell array
M_cell = num2cell(M_numeric);

% add header row
M_out = [header; M_cell];

% Export as a single table with variable names
T_export = array2table(M_out, ...
    'VariableNames', {'Component_group', ...
    'Detection_rate (%)', ...
    'Component_base_peak', ...
    '¹tᵣ Median', ...
    '²tᵣ Median', ...
    ['std in QC (n:',num2str(sum(idx_QC)),')'] ...
    ,['mean in QC (n:',num2str(sum(idx_QC)),')'] ...
    ,['RSD in QC (n:',num2str(sum(idx_QC)),')'], ...
    'blank_mean_int', ...
    'sample_min_int', ...
    'sample_max_int', ...
    matFile_list.name});
T_export(1:2,1:20)

%write to excell 
writetable(T_export,fullfile(Options.dir.results,'peak_table_final.xlsx'),'Sheet','peak_volume')