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
% mz_vec = [274.2017 ]%+ [0:2] * (12+1.007825*2)]
% mz_vec = metadata{:,"QuantificationIon_LC_"};
% compound_name = metadata{:,"Analyte"};
% include = false(length(mz_vec),1);
% compound_info = cell(length(mz_vec),8);
col = lines(length(mz_vec));
% for S = 1:length(mz_vec)
%
% mzTarget = mz_vec(S);
%
 figure(2)
feature_list = 13171;%clusters_05(6).names% find((abs(mzTarget-peak_table{:,3})./peak_table{:,3}).*1e6 < Options.ppm_dev);
% feature_list = clusters_05(14).names
% clf
S = 1
sample_vec_pca = 1:34
if ~isempty(feature_list)
    for n = 1:length(feature_list)
        %elution profile
        % clf
        subplot(1,2,1)
        [X, Y, Z_refolded, mass_spec, mz, p_tmp] = plotClusters_from_disk( ...
            feature_list(n), feature_table, clusters_all_samples_before_duplicate_removal,mzroi_aug, ...
            load_dir, Options, n_samples, col(S,:),sample_vec_pca',true);
        % mz_title = full(mz{1}(max(mass_spec{1},[],'all') == max(mass_spec{1},[],2)));
        %  p(group_idx(sample_vec_pca)) = p_tmp(1);
        % valid(group_idx(sample_vec_pca)) = true;
        view([-30 -60 20])
        xL = xlim;
        yL = ylim;
        zL = zlim;
        hold on
        if max(Z_refolded{1},[],'all')>2000

            % Mass spectrum
            subplot(1,2,2)
            % clf
            [val,a] = max(sum(mass_spec{1},2));
            stem(mz{1},sum(mass_spec{1},2)./val*100,'Color',col(S,:))
            text(mz{1},sum(mass_spec{1},2)./val*100,num2str(mz{1}),'rotation',45)
            title(round(mz{1}(a),4));
            % hold on
            xlabel('m/z');
            ylabel('Intensity (%)');
            % linkaxes
            % pause

            [r,c] = find(squeeze(sum(Z_refolded{1},[3:4])) == max(squeeze(sum(Z_refolded{1},[3:4])),[],'all'));
            % if all(abs([PesticeMix.Rt_1D_min(S),round(PesticeMix.Rt_2D_sec(S),2)]-[ X{1}(c,r),Y{1}(c,r)])<[1.35,1])
            %   subplot(1,2,1)
            [~,a] = max(sum(mass_spec{1},2));
            ppm_dev = (mz{1}(a)-mzTarget)/mzTarget*1e6;
            compound_info(S,:) = {feature_list(n),compound_name{S},'',mzTarget,full(mz{1}(a)),ppm_dev,X{1}(1,r),round(Y{1}(c,1),2)};
            %       whos Name
            %       Name = Name;
            %       pause
            %       scatter3(PesticeMix.Rt_1D_min(S),round(PesticeMix.Rt_2D_sec(S),2),max(Z_refolded{1},[],'all'),'filled')
            %       sgtitle([PesticeMix.Rt_1D_min(S),round(PesticeMix.Rt_2D_sec(S),2)])
            % % pause
            % end
            sgtitle([num2str(S),': ',metadata{S,"Analyte"}{:}])
            include(S) = true;
            potential_suspect{feature_list(n)}  = find(abs(mz{1}(a)-suspect_list{:,"expected_mz"})./mz{1}(a).*1e6<5);
            % 'stop'
            clc
        end
    end
    % 'bingo'

end

% end
%%
include2 = include;
include2([11,12,13,16,20,14,28,30,31,33,34,36,39,40,41,46]) = false
% compound_info([11,12,13,16,20,14,28,30,31,33,34,36,39,40,41],:)
%%
clc
compound_info_cleaned = compound_info(include2,:)