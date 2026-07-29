
%% --- Load Metadata ---
load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\curve_resolution_final_meta_data_Ver2.mat")
load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\curve_resolution_final_251110_Ver2.mat","feature_table")
load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\num_fragments_in_cluster.mat")
SampleList = readtable("H:\PhD\Experiments\LargeSampleSet_Vandalf\SampleList\SampleList_Combined.xlsx");
fileList_tmp = matFile_list;

% Keep only matching samples
SampleList = SampleList(contains(SampleList.SampleName_1, extractBefore({fileList_tmp.name},10)),:);

% Remove specific sample if needed
SampleList = SampleList(~contains(SampleList.SampleName_1,{'270724_38.d'}),:);

% Define groups
vec_QC     = find(ismember(SampleList.WWTP,{'QC1','QC2','QC3'}));
vec_ValQC  = find(ismember(SampleList.WWTP,{'ValQC1','ValQC2'}));
vec_blank  = find(ismember(SampleList.WWTP,{'Blank','InstrumentBlank'}));
vec_Source = find(ismember(SampleList.WWTP,{'SA','SB'}));

% Exclude groups from fileList
exclude_idx = [vec_blank; vec_QC; vec_Source; vec_ValQC];
fileList_tmp(exclude_idx) = [];

% Indices for PCA samples
sample_vec_pca = setdiff(1:size(SampleList,1), exclude_idx);

WWTP_type = unique(SampleList.WWTP);

%% --- Build Feature Matrix ---
% X = accumarray([feature_groups_all.sample; clusters_all_samples']', ...
%                [feature_groups_all.volume_summed]');
X = readtable("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\peak_table_final.xlsx","Sheet","peak_volume",'Range','L3:CK18005');
X = X{:,:}';
% Keep only features detected in >30% of samples
% include_feature = find(sum(X>0) > size(X,1)*0.3);
% X = X(:,include_feature);

% Split groups
X_blank  = X(vec_blank,:);
X_qc     = X(vec_QC,:);
X_val_qc = X(vec_ValQC,:) ;
% Set negatives to zero
X(X<1000) = 0;

% Remove excluded groups
X(exclude_idx,:) = [];
%%
% QC reproducibility filter
X_qc_std = std(X_qc,'omitnan');
rel_peaks_QC    = (X_qc_std ./ mean(X_qc,'omitnan') * 100 < 30)';  % CV <30%
rel_peaks_Blank = sum((X./mean(X_blank,'omitnan')) > 3)' > size(X,1)*0.3;
ind_keep = rel_peaks_QC & rel_peaks_Blank & num_mz' > 1;

sum(rel_peaks_QC)
sum(rel_peaks_Blank)
sum(num_mz' > 1)

% Scale by QC std
valid_idx = ~isnan(X_qc_std) & ~isinf(X_qc_std) & ind_keep';
X_scaled  = X(:,valid_idx) ./ X_qc_std(valid_idx);
% X_scaled  = (X-mean(X)) ./ std(X);%XX_qc_std(valid_idx);

% X_scaled(isnan(X_scaled) | isinf(X_scaled)) = 0;

X_qc_scaled        = X_qc(:,valid_idx) ;%./ X_qc_std(valid_idx);
X_val_qc_scaled    = X_val_qc(:,valid_idx) ;%./ X_qc_std(valid_idx);
final_features = find(valid_idx);
% X_pca = X(:,valid_idx);
%% --- PCA ---
[X_pca,mu,sigma] = zscore(X(:,valid_idx));
% mu = mean(X_pca);
% std_X = std(X_pca);
% % (X_pca-mu)./std(X)

[coeff, score, latent, tsquared, explained,~] = pca(X_pca, 'Centered', 'off');
% final_features = include_feature(valid_idx);
[~,a] = max(abs(coeff));
vec_var  = [1:size(coeff,1)];
% final_features = final_features(~ismember(vec_var,a(1:2)));
% vec_var = vec_var(~ismember(vec_var,a(1:2)));

% [coeff, score, latent, tsquared, explained, mu] = pca(X_scaled(:,vec_var));

%% --- PCA Scores Plot ---

figure(1)
clf
subplot(2,2,1)
pc_select = [1,2];
col = lines(length(WWTP_type));
col(end-5,[2,3]) = [0.5,1];
col(end-2,:) = [1,0.5,0.4];
col(2,:) = [0.5660    0.4430    0.7450];
hold on
legend_entries = {};
legend_handles = [];
SampleList_excluded = SampleList;
SampleList_excluded(exclude_idx,:) = [];
for g = 1:length(WWTP_type)
    ind = strcmp(SampleList.WWTP, WWTP_type{g});
    ind(exclude_idx) = [];  % remove excluded samples
    groupScores = score(ind,:);

    if isempty(groupScores), continue; end

    % Scatter points
    h = scatter(groupScores(:,pc_select(1)), groupScores(:,pc_select(2)), 50, col(g,:), 'filled');
labels = string(SampleList_excluded{ind,'StartDate'});
labels = extractBefore(labels, 7);   % first 6 characters

text(groupScores(:,pc_select(1)), groupScores(:,pc_select(2)), labels)
   
    legend_handles(end+1) = h; %#ok<AGROW>
    legend_entries{end+1} = WWTP_type{g}; %#ok<AGROW>

    % Convex hull for groups (except QC)
    if size(groupScores,1) > 2 && ~ismember(WWTP_type{g},{'QC1','QC2','QC3'})
        K = convhull(groupScores(:,pc_select(1)), groupScores(:,pc_select(2)));
        patch(groupScores(K,pc_select(1)), groupScores(K,pc_select(2)), col(g,:), ...
            'FaceAlpha',0.15, 'EdgeColor',col(g,:), 'LineWidth',1.5);
    end
end

% Project QC samples onto PCA space
X_qc_centered = (X_qc(:,valid_idx)  - mu)./sigma;
qc_scores = X_qc_centered * coeff;
legend_handles(end+1) = scatter(qc_scores(:,pc_select(1)), qc_scores(:,pc_select(2)), 50, 'r', 'filled');
legend_entries =[legend_entries,{'QC'}];
legend(legend_handles, legend_entries, 'Location','northeast')

% Project QC samples onto PCA space
X_val_qc_centered = (X_val_qc_scaled(:,vec_var) - mu)./sigma;
X_val_qc_scores = X_val_qc_centered * coeff;

legend_handles(end+1) = scatter(X_val_qc_scores(:,pc_select(1)), X_val_qc_scores(:,pc_select(2)), 50, 'b', 'filled');
legend_entries =[legend_entries,{'Val QC'}];
legend(legend_handles, legend_entries, 'Location','northeast')


axis equal
xlabel(['PC',num2str(pc_select(1)),' (', num2str(round(explained(pc_select(1)),1)), '%)'])
ylabel(['PC',num2str(pc_select(2)),' (', num2str(round(explained(pc_select(2)),1)), '%)'])

text(0.05,0.93,[figLetter(1),') PCA Scores'],'units','normalized','FontSize',25)
% set(gca,'XDir','reverse')
run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')
% --- PCA Loadings Plot ---
subplot(2,2,3)
PCA_matching_suspects
% match_IDs
% scatter(coeff(:,pc_select(1)), coeff(:,pc_select(2)),'b','filled')
% t = text(coeff(:,1), coeff(:,2), num2str(final_features'));
xlabel(['PC',num2str(pc_select(1)),' (', num2str(round(explained(pc_select(1)),1)), '%)'])
ylabel(['PC',num2str(pc_select(2)),' (', num2str(round(explained(pc_select(2)),1)), '%)'])
text(0.05,0.93,[figLetter(2),') PCA Loadings'],'units','normalized','FontSize',25)
run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')
% set(gca,'XDir','reverse')
%% --- Example Feature Visualization ---
load_dir = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\component_cluster\';
load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob\nnmf_251119_ver3\curve_resolution_final_251110_Ver2.mat","clusters_all_samples_before_duplicate_removal")
% figure(1)
clc
% clf
feature_list =[980,2836,2388,2932,4057];[985,1029,1090,720,1352];[1760,1553,1780,1809,1753];[2316,2122,2289,2178,1869];
% for  nn = 1:length(feature_list)
%     subplot(2,2,3)
% %    ;
%     ind_f = final_features== feature_list(nn); % ind_f = include_feature(valid_idx)== feature_list(nn);
%    if sum(ind_f)>0
%         t(ind_f).Color = 'r';
%    end
% end
% %
for k = 1:78
    group_idx(k) = find(strcmp(WWTP_type, SampleList.WWTP{k}));
end

L_compound_name = {'cotinine','4-amino-5-hydroxy-3-oxotetradecanoic acid','lamotrigine','unidentified ','clopidogrel'}
for nn = 1:length(feature_list)
    fn = feature_list(nn);
    subplot(5,2,nn*2)%ceil(sqrt(length(feature_list))), ceil(sqrt(length(feature_list)))-1, nn)
    hold on
    valid = false(length(WWTP_type),1);
    clear p
    % for k = sample_vec_pca
    %     a = find(clusters_all_samples == fn & [feature_groups_all.sample]' == k,1,'first');
    %     if isempty(a), continue; end
    %
    %     group_idx = find(strcmp(WWTP_type, SampleList.WWTP{k}));
    %     p_tmp = plotFeatureGroup3D(feature_groups_all, Options, a, 'color', col(group_idx,:));
    %     zlabel('Intensity (a.u)')
    %     valid(group_idx) = true;
    %     p(group_idx) = p_tmp(1);
    %     mz_title = full(feature_groups_all(a).basePeak);
    %   
    % end
    
    [~, Y, Z_refolded, mass_spec, mz, p_tmp] = plotClusters_from_disk( ...
        feature_list(nn), feature_table, clusters_all_samples_before_duplicate_removal,mzroi_aug, ...
        load_dir, Options, 78, col(group_idx,:),sample_vec_pca,false);
    mz_title = full(mz{1}(max(mass_spec{1},[],'all') == max(mass_spec{1},[],2)));
    % p(group_idx(sample_vec_pca)) = p_tmp(1);
    % valid(group_idx(sample_vec_pca)) = true;
    view([-30 -60 20])
    xL = xlim;
    yL = ylim;
    zL = zlim;
    % legend(p(valid),WWTP_type(valid), 'Location','northeastoutside')
    
    % text(xL(1),yL(2),zL(2),[sprintf([figLetter(nn+2),') m/z: %.4f | '], mz_title), L_compound_name{nn}],'FontSize',10,'BackgroundColor','W','color','k')
    text(0.02,0.9,1,[sprintf([figLetter(nn+2),') m/z: %.4f | '], mz_title), L_compound_name{nn}],'FontSize',10,'BackgroundColor','W','color','k','Units','normalized')
   
    %    run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')
    set(gca,'FontSize',8,'LineWidth',2)
    box on
end

%%
fig = gcf;
fig.Position =   [1        1        1278         946];
exportgraphics(fig,'H:\PhD\Articles\8_LCxLC_Vandalf_BigDataSet\Part_2_Data\Figures\PCA_NTS.tif','Resolution',1000)

