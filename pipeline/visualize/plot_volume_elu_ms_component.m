function volume = plot_volume_elu_ms_component( ...
    num_clusters_before_curve_resolution,clusters, ...                         % <-- directly provided
    feature_table, rt_offset, Options, ...
    clusters_all_samples_before_duplicate_removal, ...
    mzroi_aug, load_dir, n_samples, ...
    sample_names, compound_name)

clf

%% ================== GROUPING ==================
[clusters_unique,~,group_idx] = unique(clusters);

vec = ismember(feature_table(:,2), clusters_unique);

% Recompute grouping based on selected clusters
[clusters,~,group_idx] = unique(feature_table(vec,2));

volume = accumarray([feature_table(vec,5), group_idx], ...
                    feature_table(vec,6), [n_samples, max(group_idx)]);

Rt_tmp = [ ...
    accumarray(group_idx, feature_table(vec,3), [], @median), ...
    accumarray(group_idx, feature_table(vec,4), [], @median)];

Rt_tmp = round(Rt_tmp);

colors = rand(max(group_idx),3);

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
        1:n_samples, false);

    % ---- Find global max in Z ----
    Z = squeeze(max(Z_refolded{1},[],3:4));   % extract matrix

    [r,c] = find(Z == max(Z(:)));
    % [row, col] = ind2sub(size(Z), idx_max);
    %
    % ---- Convert to axis coordinates ----
    x_peak = X_coord{1}(1, r);
    y_peak = Y_coord{1}(2);

    % ---- Add label ----
    text(x_peak, y_peak, max(Z(:)), num2str(n_cluster), ...
        'Color', colors(k,:), ...
        'FontWeight','bold', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom')
    % volume(:,k) =  squeeze(sum(Z_refolded{1},1:3));
    %% ---- Spectrum ----

    subplot(2,2,4)
    y_spec = sum(mass_spec{n_cluster},2);

    stem(mz{n_cluster}, y_spec, ...
        'Marker','none','Color',colors(k,:))
    ind = y_spec>0.01*max(y_spec);
    text(mz{n_cluster}(ind), y_spec(ind), ...
        num2str(mz{n_cluster}(ind)), ...
        'Rotation',45,'Color',colors(k,:))

    xlabel('m/z')
    hold on
end

%% ================== AXIS FORMATTING ==================
subplot(2,2,2)
xt = xticks;
% xticklabels(string(xt + rt_offset))

subplot(2,2,4)
legend(string(clusters))
% %% ================== BAR PLOT ==================
% subplot(1,2,1); hold on
% 
% b = gobjects(max(group_idx),1);
% 
% for k = 1:max(group_idx)
%     b(k) = bar(volume(:,k), 'FaceColor', colors(k,:));
% end
% 
% xticks(1:length(sample_names))
% xticklabels(sample_names)
% ylabel('Volume')
% legend(b, string(clusters))

%% ---- Axis formatting ----
subplot(2,2,2)
xt = xticks;
% xticklabels(string(xt + rt_offset))

subplot(2,2,4)
legend(string(clusters))

%% ---- Save name ----
% name = regexprep(compound_name{j}, '[<>:"/\\|?*]', '_');

%% ================== TITLE ==================
if exist('compound_name','var') && ~isempty(compound_name)
    sgtitle(compound_name)
end

end