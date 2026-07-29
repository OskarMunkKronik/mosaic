function plot_elution_profile_and_mass_spectrum( ...
    feature_list, feature_table, clusters_all_samples_before_duplicate_removal, mzroi_aug, ...
    load_dir, Options, n_samples, col,  sample_vec_pca)

% plotFeatureClusters
% Plots elution profile and mass spectrum for a list of features
%
% Inputs:
%   feature_list      - vector of feature indices
%   feature_table     - feature table
%   clusters_all_samples - cluster structure
%   mzroi_aug         - m/z ROI
%   load_dir          - data directory
%   Options           - options struct
%   n_samples         - number of samples
%   col               - color matrix
%   S                 - sample index
%   sample_vec_pca    - PCA sample vector

    if isempty(feature_list)
        return;
    end

    for n = 1:length(feature_list)

        %% Elution profile
        subplot(length(feature_list),2,n*2-1)
        [X,Y,Z_refolded, mass_spec, mz, ~] = plotClusters_from_disk( ...
            feature_list(n), feature_table, clusters_all_samples_before_duplicate_removal,mzroi_aug, ...
            load_dir, Options, n_samples, col,sample_vec_pca,true);
        [r,c] = find(squeeze(sum(Z_refolded{1},[3:4])) == max(squeeze(sum(Z_refolded{1},[3:4])),[],'all'));
            title(X{1}(1,r),round(Y{1}(c,1),2))
        view([-30 -60 20])
        hold on

        %% Mass spectrum
            subplot(length(feature_list),2,n*2)
            for k = 1:n_samples
        spectrum_sum = sum(mass_spec{1}(:,k), 2);
        [val, a] = max(spectrum_sum);

        intensity_norm = spectrum_sum ;%./ val * 100;

        stem(mz{1}, intensity_norm, 'Color', col(k,:))
        hold on
        end 
        text(mz{1}, max(mass_spec{1},[],2), num2str(mz{1}), 'Rotation', 45)
            
        title(round(mz{1}(a), 4))
        xlabel('m/z')
        ylabel('Intensity (%)')

    end
end
