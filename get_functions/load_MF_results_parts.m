function [clusters_all_samples, feature_groups_all] = load_MF_results_parts(savePath,prefix)
%LOAD_MF_RESULTS_PARTS Load MF results saved in multiple .mat files.
%
%   [clusters_all_samples, feature_groups_all] = ...
%       load_MF_results_parts(savePath)
%
% INPUT:
%   savePath - folder containing:
%       MF_results_part_001.mat
%       MF_results_part_002.mat
%       MF_results_part_003.mat
%       ...
%
% OUTPUT:
%   clusters_all_samples
%   feature_groups_all
%
% The function reconstructs the original cell arrays by concatenating
% the parts in numerical order.

    % -------------------------------------------------------------
    % Find all parts
    % -------------------------------------------------------------

    files = dir(fullfile(savePath, [prefix,'*.mat']));

    if isempty(files)
        error('No MF_results_part_*.mat files found in:\n%s', savePath);
    end

    % -------------------------------------------------------------
    % Sort files by part number
    % -------------------------------------------------------------

    filenames = {files.name};

    partNumbers = zeros(numel(filenames), 1);

    for i = 1:numel(filenames)

        token = regexp( ...
            filenames{i}, ...
            [prefix,'(\d+)\.mat'], ...
            'tokens', ...
            'once');

        if isempty(token)
            error('Could not determine part number from: %s', ...
                filenames{i});
        end

        partNumbers(i) = str2double(token{1});
    end

    [~, order] = sort(partNumbers);

    files = files(order);

    % -------------------------------------------------------------
    % Initialize
    % -------------------------------------------------------------

    clusters_all_samples = [];
    feature_groups_all = struct([]);

    fprintf('Found %d MF result files.\n', numel(files));

    % -------------------------------------------------------------
    % Load each part
    % -------------------------------------------------------------

    for p = 1:numel(files)

        filename = fullfile(files(p).folder, files(p).name);

        fprintf('Loading part %d/%d: %s\n', ...
            p, numel(files), files(p).name);

        S = load(filename, ...
            'clusters_part', ...
            'feature_groups_part');

        % ---------------------------------------------------------
        % Append to output
        % ---------------------------------------------------------

        clusters_all_samples = cat(1,clusters_all_samples(:), S.clusters_part(:));

        if isempty(feature_groups_all)
            feature_groups_all = S.feature_groups_part(:);
        else

        feature_groups_all = cat(1, feature_groups_all,S.feature_groups_part(:));
        end 
      
        clear S

    end

    % -------------------------------------------------------------
    % Check result
    % -------------------------------------------------------------

    fprintf('\nSuccessfully loaded all parts.\n');
    fprintf('Total samples: %d\n', ...
        numel(clusters_all_samples));

end