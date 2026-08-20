function save_MF_results_parts(clusters_all_samples, feature_groups_all, savePath,prefix, chunkSize)
%SAVE_MF_RESULTS_PARTS Save large MF results into multiple .mat files.
%
% Automatically reduces chunk size if an out-of-memory error occurs.
%
% INPUTS:
%   clusters_all_samples - cell array, one entry per sample
%   feature_groups_all   - cell array, one entry per sample
%   savePath             - folder where files should be saved
%   chunkSize             - initial number of samples per file

    % -------------------------------------------------------------
    % Input checks
    % -------------------------------------------------------------

    if ~isempty(clusters_all_samples) && ...
        numel(clusters_all_samples) ~= numel(feature_groups_all)
    error(['clusters_all_samples and feature_groups_all ' ...
           'must have the same number of samples.']);
end

    if nargin < 4 || isempty(chunkSize)
        chunkSize = 10;
    end

    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end

    nSamples = numel(clusters_all_samples);
    chunkSize = ceil(nSamples/chunkSize);
    % -------------------------------------------------------------
    % Try progressively smaller chunk sizes
    % -------------------------------------------------------------

    while chunkSize >= 1

        fprintf('\nTrying chunk size = %d samples/file\n', chunkSize);

        try

            % IMPORTANT:
            % Number of files depends on total samples,
            % not on chunkSize alone.
            nParts = ceil(nSamples / chunkSize);

            for p = 1:nParts

                idxStart = (p-1)*chunkSize + 1;
                idxEnd   = min(p*chunkSize, nSamples);

                idx = idxStart:idxEnd;

                fprintf('  Saving part %d/%d (samples %d-%d)...\n', ...
                    p, nParts, idxStart, idxEnd);

                % Extract current chunk
                % if 
                clusters_part = clusters_all_samples(idx);
                feature_groups_part = feature_groups_all(idx);

                % Filename
                filename = fullfile( ...
                    savePath, ...
                    sprintf('%s_%03d.mat', prefix, p));
                % Save
                save(filename, ...
                    'clusters_part', ...
                    'feature_groups_part');

                % Clear temporary variables
                clear clusters_part feature_groups_part idx

            end

            % -----------------------------------------------------
            % Everything succeeded
            % -----------------------------------------------------

            fprintf('\nSuccessfully saved all data.\n');
            fprintf('Chunk size: %d samples/file\n', chunkSize);
            fprintf('Number of files: %d\n', nParts);

            return

        catch ME

            % Clear temporary variables
            clear clusters_part feature_groups_part idx

            % Check if this was an out-of-memory error
            isOutOfMemory = ...
                contains(ME.message, 'Out of memory', ...
                         'IgnoreCase', true) || ...
                contains(ME.identifier, 'MATLAB:nomem', ...
                         'IgnoreCase', true);

            if isOutOfMemory

                fprintf('\nOUT OF MEMORY with chunk size %d.\n', ...
                    chunkSize);

                % Delete potentially incomplete file
                if exist('filename', 'var') && isfile(filename)
                    delete(filename);
                end

                % Reduce chunk size
                newChunkSize = max(1, floor(chunkSize / 2));

                if newChunkSize == chunkSize
                    error(['Unable to save even one sample without ' ...
                           'running out of memory.']);
                end

                chunkSize = newChunkSize;

                fprintf('Retrying with chunk size = %d...\n\n', ...
                    chunkSize);

            else
                % Not a memory error -- don't hide the real problem
                rethrow(ME);
            end
        end
    end

    error('Could not save the data.');

end