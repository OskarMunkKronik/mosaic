function matches = search_plot_volume_elu_ms( ...
    num_clusters_before_curve_resolution,mz_in_mass_spectra, mzTarget, ppm_dev, ...
    feature_table, rt_offset, Options, ...
    clusters_all_samples_before_duplicate_removal, ...
    mzroi_aug, load_dir, n_samples, ...
    sample_names, compound_name)

clf

%% ================== MATCHING ==================
nFeatures = numel(mz_in_mass_spectra);
nTargets  = size(mzTarget,1);
% Options.coords.X = Options.coords.X + rt_offset;
        

matches = false(nFeatures, nTargets);

for i = 1:nFeatures
    mz_spec = mz_in_mass_spectra{i};

    for j = 1:nTargets
        ppm_error = abs(mz_spec - mzTarget(j,:)) ./ mzTarget(j,:) * 1e6;
        matches(i,j) = any(sum(ppm_error < ppm_dev,2) > 0);
    end
end

% disp(sum(matches))

%% ================== RT OFFSET ==================
% rt_offset = round(median(cellfun(@(rt) rt(1),Rt_aug))/60,3);

%% ================== LOOP TARGETS ==================
for j = 1:nTargets
    
    vec = find(matches(:,j));
    if isempty(vec)
        continue
    end
    
    clf
    
    %% ---- Grouping ----
    [clusters,~,group_idx] = unique(feature_table(vec,2));
    
    volume = accumarray([feature_table(vec,5), group_idx], ...
                        feature_table(vec,6), [n_samples, max(group_idx)]);
    
    Rt_tmp = [ ...
        accumarray(group_idx, feature_table(vec,3), [], @median), ...
        accumarray(group_idx, feature_table(vec,4), [], @median)];
    
    Rt_tmp = round(Rt_tmp);

    colors = rand(max(group_idx),3);

 
    %% ================== CHROMATOGRAM + SPECTRUM ==================
    clear mass_spec mz
    
    for k = 1:numel(clusters)
        n_cluster = clusters(k);
        
        %% ---- Chromatogram ----
        subplot(2,2,2)
        [X_coord,Y_coord,Z_refolded,mass_spec(n_cluster),mz(n_cluster),~] = ...
            plotClusters_from_disk( ...
            num_clusters_before_curve_resolution,n_cluster, feature_table, ...
            clusters_all_samples_before_duplicate_removal', ...
            mzroi_aug, load_dir, Options, n_samples, ...
            repmat(colors(k,:),[n_samples,1]), ...
            1:n_samples, true);

        % ---- Find global max in Z ----
        Z = squeeze(max(Z_refolded{1},[],3:4));   % extract matrix
        
        [r,c] = find(Z == max(Z(:)));
        % [row, col] = ind2sub(size(Z), idx_max);
        % 
        % ---- Convert to axis coordinates ----
        x_peak = X_coord{1}(1, r);
        y_peak = Y_coord{1}(c);
        
        % ---- Add label ----
        if  max(Z(:))> 0 
        text(x_peak, y_peak, max(Z(:)), num2str(n_cluster), ...
            'Color', colors(k,:), ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom')
        end 

        %% ---- Spectrum ----
        subplot(2,2,4)
        y_spec = sum(mass_spec{n_cluster},2);
        
        stem(mz{n_cluster}, y_spec, ...
            'Marker','none','Color',colors(k,:))
        
        text(mz{n_cluster}, y_spec, ...
            num2str(mz{n_cluster}), ...
            'Rotation',45,'Color',colors(k,:))
        
        xlabel('m/z')
        hold on

        volume(:,k) =  squeeze(sum(Z_refolded{1},1:3));
    end
       %% ================== BAR PLOT ==================
    subplot(1,2,1); hold on
    
    b = gobjects(max(group_idx),1);
    
    for k = 1:max(group_idx)
        b(k) = bar(volume(:,k), 'FaceColor', colors(k,:));
    end
    
    xticks(1:length(sample_names))
    xticklabels(sample_names)
    ylabel('Volume')
    legend(b, string(clusters))

    %% ---- Axis formatting ----
    subplot(2,2,2)
    xt = xticks;
    % xticklabels(string(xt + rt_offset))

    subplot(2,2,4)
    legend(string(clusters))

    %% ---- Save name ----
    name = regexprep(compound_name{j}, '[<>:"/\\|?*]', '_');

    % Optional saving:
    % savefig(gcf, fullfile('path_to_save',[name '.fig']))
end

end