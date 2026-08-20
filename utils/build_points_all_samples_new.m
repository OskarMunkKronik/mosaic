function points = build_points_all_samples_new(fileList, Options, mzroi_aug)
%BUILD_POINTS_ALL_SAMPLES
% Build feature matrix/cell array from clustering and peak detection results.
%
% Each ROW = one feature
% Each COLUMN = one variable
%
% Columns:
%   1  Rt
%   2  elution1
%   3  elution2
%   4  feature_elu_2D
%   5  massSpec
%   6  volume_summed
%   7  volume_basePeak
%   8  basePeak
%   9  sample
%   10 sample_feature_ID
%   11 sample_cluster_ID
%   12 IS_flag
%
% Matrix-compatible variables are stored directly.
% Variable-size variables are stored in cells.

    % -------------------------------------------------------------
    % Column definitions
    % -------------------------------------------------------------
    col.Rt                = 1;
    col.elution1          = 2;
    col.elution2          = 3;
    col.feature_elu_2D    = 4;
    col.massSpec          = 5;
    col.volume_summed     = 6;
    col.volume_basePeak   = 7;
    col.basePeak          = 8;
    col.sample             = 9;
    col.sample_feature_ID  = 10;
    col.sample_cluster_ID  = 11;
    col.IS_flag            = 12;

    nCols = 12;

    % -------------------------------------------------------------
    % Preallocate as cell array
    % -------------------------------------------------------------
    points = cell(0, nCols);

    % Header describing columns
    points_header = { ...
        'Rt', ...
        'elution1', ...
        'elution2', ...
        'feature_elu_2D', ...
        'massSpec', ...
        'volume_summed', ...
        'volume_basePeak', ...
        'basePeak', ...
        'sample', ...
        'sample_feature_ID', ...
        'sample_cluster_ID', ...
        'IS_flag'};

    % -------------------------------------------------------------
    % Process samples
    % -------------------------------------------------------------
    row = 0;

    for k = 1:length(fileList)

        % Find MF file
        files = dir(fullfile( ...
            Options.Paths.save2mat, ...
            'MF', ...
            [fileList(k).name '*.mat']));

        if isempty(files)
            error('No matching file found for %s*.mat', ...
                fileList(k).name);
        end

        % Load file
        load(fullfile(files(1).folder, files(1).name));

        % ---------------------------------------------------------
        % Find valid features
        % ---------------------------------------------------------
        valid = allPeaks(:,6) > 0;
        feature_idx = find(valid);

        % ---------------------------------------------------------
        % Process each feature
        % ---------------------------------------------------------
        for ii = 1:length(feature_idx)

            n = feature_idx(ii);

            row = row + 1;

            % Cluster ID
            cluster_id = clusters(n);

            % Features belonging to this cluster
            cluster_mask = (clusters == cluster_id);

            % -----------------------------------------------------
            % Basic feature information
            % -----------------------------------------------------

            points{row, col.Rt} = allPeaks(n,[2 3])';

            points{row, col.elution1} = sum_elution_profile{n}{1};

            points{row, col.elution2} = sum_elution_profile{n}{2};

            points{row, col.feature_elu_2D} = Z_feature(n);

            points{row, col.massSpec} = mass_spectra(:,cluster_id);

            % -----------------------------------------------------
            % Quantities
            % -----------------------------------------------------

            points{row, col.volume_summed} = ...
                sum(allPeaks(cluster_mask,5));

            [points{row, col.volume_basePeak}, a] = ...
                max(mass_spectra(:,cluster_id));

            points{row, col.basePeak} = mzroi_aug(a);

            % -----------------------------------------------------
            % Metadata
            % -----------------------------------------------------

            points{row, col.sample} = k;

            points{row, col.sample_feature_ID} = ...
                allPeaks(n,6);

            points{row, col.sample_cluster_ID} = ...
                cluster_id;

            % -----------------------------------------------------
            % IS flag
            % -----------------------------------------------------

            if size(allPeaks,2) > 6
                points{row, col.IS_flag} = allPeaks(n,7);
            else
                points{row, col.IS_flag} = NaN;
            end

        end
    end

    % -------------------------------------------------------------
    % Store header as first row if desired
    % -------------------------------------------------------------
    % points = [points_header; points];

end