function [X, Y, Z_refolded, mass_spec, mz_vals_cluster,p] = plotClusters_from_disk( ...
    num_clusters_before_curve_resolution,n_cluster_vec, feature_table, clusters_all_samples_before_duplicate_removal,mzroi_aug, ...
    load_dir, Options, n_samples, col,sample_vec,faster)

% plotClusters
% Loops through cluster IDs, loads results files, refolds components,
% generates 3D plots and mass spectra, and returns key variables.
%
% INPUTS:
%   n_cluster_vec   – vector of cluster numbers to process
%   feature_table   – feature table, column 2 holding cluster IDs
%   clusters_all_samples_before_duplicate_removal – cluster index vector
%   load_dir        – folder containing 'results_XXX.mat' files
%   Options         – struct containing coords.X and coords.Y
%   n_samples       – number of samples for refolding
%   col             – Nx3 color matrix for plotting
%
% OUTPUTS:
%   X, Y            – cell arrays of meshgrid coordinates
%   Z_refolded      – refolded component data
%   mass_spec       – mass spectra for each cluster
%   mz_vals_cluster – mz values for each cluster
%
% ---------------------------------------------------------------------

n_clusters = length(n_cluster_vec);

% Preallocate outputs
X = cell(n_clusters, 1);
Y = cell(n_clusters, 1);
Z_refolded = cell(n_clusters, 1);
mass_spec = cell(n_clusters, 1);
mz_vals_cluster = cell(n_clusters, 1);

% figure;   % one figure containing both subplots
% cd(load_dir)
% numFiles  = numel(dir('*.mat'));
% n_digits_clusters = length(num2str(numFiles));
n_digits_clusters = ceil(log10(num_clusters_before_curve_resolution )+ 1);

for n = 1:n_clusters
    
    % -----------------------------------------
    % Select samples belonging to this cluster
    % -----------------------------------------
    ind_cluster = feature_table(:,2) == n_cluster_vec(n);

    % Find cluster group with most counts
    clusters_all_samples_before_duplicate_removal = clusters_all_samples_before_duplicate_removal(:);
    [gc, grps] = groupcounts(clusters_all_samples_before_duplicate_removal(ind_cluster));
    n_cluster_old = grps(1);

    % -----------------------------------------
    % Load results file
    % -----------------------------------------

    % Zero-padded cluster index
    cluster_str = sprintf(['%0', num2str(n_digits_clusters), 'd'], n_cluster_old);

    % Filename: results_XXX.mat
    filename = sprintf("results_%s.mat", cluster_str);

    % Load data
    results = load(fullfile(load_dir, filename));
  
    % -----------------------------------------
    % Generate X,Y meshgrid
    % -----------------------------------------
    [X{n}, Y{n}] = meshgrid( ...
        Options.coords.X(1, results.x_range), ...
        Options.coords.Y(results.y_range, 1));

    % -----------------------------------------
    % Refold component images
    % -----------------------------------------
    % Z_refolded{n} = refold_Z_blocks( ...
    %     results.W, results.H, ...
    %     results.x_range, results.y_range, results.mz_range, ...
    %     n_samples, []);

    % Obtain mz-values
    mz_vals_cluster{n} = mzroi_aug(results.mz_range);

    % -----------------------------------------
    % PLOT: Spatial slices (subplot 1)
    % -----------------------------------------
    % subplot(2,1,1);
    hold on;

    mass_spec{n} = zeros(length(mz_vals_cluster{n}), size(results.Z_refolded,4));
    if faster

        Zk = squeeze(max(results.Z_refolded,[],3:4));

        p = plot3(X{n}, Y{n}, Zk, 'Color', col(1,:),'LineWidth',1.5);

        % p = surf(X{n}, Y{n}, Zk','EdgeColor','interp','FaceColor','interp');
        % mass spectrum
        mass_spec{n} = zeros(size(results.Z_refolded,[3,4]));
        mass_spec{n}(:,:) = squeeze(sum(results.Z_refolded, [1 2]));
        n
        view([-10 -30 20])
    else
        for k_count = 1:length(sample_vec)%1:size(results.Z_refolded{n},4)
            k = sample_vec(k_count);
            % 3D plot of max-projected spatial component
            Zk = squeeze(max(results.Z_refolded(:,:,:,k),[],3));
            p = plot3(X{n}, Y{n}, Zk, 'Color', col(k,:),'LineWidth',1.5);
            hold on

            % mass spectrum
            mass_spec{n}(:,k) = squeeze(sum(results.Z_refolded(:,:,:,k), [1 2]));
        end
        view([-10 -30 20])
    end
    xlabel('^1t_R (min)')

    ylabel('^2t_R (sec)')
    zlabel('Intensity')
    % -----------------------------------------
    % PLOT: Spectrum (subplot 2)
    % -----------------------------------------
    % % subplot(2,1,2);
    % hold on;
    % stem(mz_vals_cluster{n}, median(mass_spec{n},2), 'Color', col(n,:));

end
Z_refolded{1} =results.Z_refolded;
end
