
%% Calculate for table 
T = accumarray([feature_table(:,2),feature_table(:,5)],feature_table(:,6));
Rt = zeros(max(feature_table(:,2)),2);
Rt(:,1) = accumarray(feature_table(:,2),Options.coords.X(1,feature_table(:,3)),[],@median);
Rt(:,2) = accumarray(feature_table(:,2),Options.coords.Y(feature_table(:,4),1),[],@median);
basepeak = accumarray(feature_table(:,2),feature_table(:,1),[],@median);
Component_group = accumarray(feature_table(:,2),feature_table(:,2),[],@median);
detection_rate_all_samples =round(sum(T>0,2)*100./size(T,2));
% idx_QC = contains(sample_id,'QC') & ~contains(sample_id,'ValQC') ;
% QC_std = std(T(:,idx_QC),[],2);
% QC_mean = mean(T(:,idx_QC),2);
% rsd_qc = QC_std./QC_mean .* 100;
% blank_mean = mean(T(:,contains(sample_id,'Blank') & ~contains(sample_id,'InstrumentBlank')),2);
% max_int = max(T(:,~contains(sample_id,{'QC','Blank'})) ,[],2);
% min_int = min(T(:,~contains(sample_id,{'QC','Blank'})) ,[],2);
% Combine into one matrix

%% Export to table 
% clc


% numeric part
M_numeric = [Component_group, detection_rate_all_samples, basepeak, Rt,T];

% header row as strings
header = [{'#', '%', 'm/z', '¹tᵣ (min)', '²tᵣ (sec)'}, {  fileList.name}];

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
      fileList.name});

%write to excell 
writetable(T_export,fullfile(Options.Paths.save2mat,'peak_table_final.xlsx'),'Sheet','peak_volume')

