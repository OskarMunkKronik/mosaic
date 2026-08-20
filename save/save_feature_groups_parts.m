
function save_feature_groups_parts(feature_groups_all, savePath, prefix, chunkSize)
%SAVE_FEATURE_GROUPS_PARTS Save feature_groups_all into multiple .mat files.
%
% Automatically reduces chunk size if an out-of-memory error occurs.
%
% INPUTS:
%   feature_groups_all - structure array, one entry per sample
%   savePath           - folder where files should be saved
%   prefix             - prefix for the output filenames
%   chunkSize          - number of samples per file
%
% Example:
%   save_feature_groups_parts(feature_groups_all, savePath, ...
%                             'MF_results', 10)

% -------------------------------------------------------------
% Input checks
% -------------------------------------------------------------

if nargin < 2 || isempty(savePath)
    error('savePath must be specified.');
end

if nargin < 3 || isempty(prefix)
    prefix = 'feature_groups';
end

if nargin < 4 || isempty(chunkSize)
    chunkSize = 10;
end

if ~isstruct(feature_groups_all)
    error('feature_groups_all must be a structure array.');
end

if ~exist(savePath, 'dir')
    mkdir(savePath);
end

nSamples = numel(feature_groups_all);
chunkSize = ceil(nSamples/chunkSize);

% Nothing to save
if nSamples == 0
    fprintf('feature_groups_all is empty. Nothing to save.\n');
    return;
end

% -------------------------------------------------------------
% Try progressively smaller chunk sizes
% -------------------------------------------------------------

while chunkSize >= 1

    fprintf('\nTrying chunk size = %d samples/file\n', chunkSize);

    try

        nParts = ceil(nSamples / chunkSize);

        for p = 1:nParts

            idxStart = (p-1)*chunkSize + 1;
            idxEnd   = min(p*chunkSize, nSamples);

            idx = idxStart:idxEnd;

            fprintf('  Saving part %d/%d (samples %d-%d)...\n', ...
                p, nParts, idxStart, idxEnd);

            % Extract current chunk
            feature_groups_part = feature_groups_all(idx);

            % Filename
            filename = fullfile( ...
                savePath, ...
                sprintf('%s_%03d.mat', prefix, p));

            % Save only the current structure chunk
            save(filename, 'feature_groups_part');

            % Clear temporary variables
            clear feature_groups_part idx

        end

        % -----------------------------------------------------
        % Everything succeeded
        % -----------------------------------------------------

        fprintf('\nSuccessfully saved all data.\n');
        fprintf('Chunk size: %d samples/file\n', chunkSize);
        fprintf('Number of files: %d\n', nParts);

        return

    catch ME

        clear feature_groups_part idx

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
