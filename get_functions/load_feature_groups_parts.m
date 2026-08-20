function feature_groups_all = load_feature_groups_parts(savePath, prefix)
%LOAD_FEATURE_GROUPS_PARTS Load feature_groups_all from multiple .mat files.
%
% INPUTS:
%   savePath - folder containing the .mat files
%   prefix   - prefix used when saving the files
%
% OUTPUT:
%   feature_groups_all - combined structure array
%
% Example:
%   feature_groups_all = load_feature_groups_parts(savePath, 'MF_results');

% -------------------------------------------------------------
% Input checks
% -------------------------------------------------------------

if nargin < 1 || isempty(savePath)
    error('savePath must be specified.');
end

if nargin < 2 || isempty(prefix)
    prefix = 'feature_groups';
end

if ~exist(savePath, 'dir')
    error('Save path does not exist: %s', savePath);
end

% -------------------------------------------------------------
% Find all matching files
% -------------------------------------------------------------

files = dir(fullfile(savePath, sprintf('%s_*.mat', prefix)));

if isempty(files)
    error('No files found with prefix "%s" in:\n%s', ...
        prefix, savePath);
end

% -------------------------------------------------------------
% Sort files by part number
% -------------------------------------------------------------

fileNames = {files.name};

partNumbers = nan(numel(fileNames), 1);

for k = 1:numel(fileNames)

    tokens = regexp(fileNames{k}, ...
        sprintf('^%s_(\\d+)\\.mat$', regexptranslate('escape', prefix)), ...
        'tokens');

    if ~isempty(tokens)
        partNumbers(k) = str2double(tokens{1}{1});
    end
end

% Remove files that don't match the expected format
valid = ~isnan(partNumbers);

files = files(valid);
partNumbers = partNumbers(valid);

[~, order] = sort(partNumbers);

files = files(order);

if isempty(files)
    error('No valid part files found.');
end

% -------------------------------------------------------------
% Load all parts
% -------------------------------------------------------------

feature_groups_all = struct([]);

for p = 1:numel(files)

    filename = fullfile(savePath, files(p).name);

    fprintf('Loading part %d/%d: %s\n', ...
        p, numel(files), files(p).name);

    data = load(filename, 'feature_groups_part');

    if ~isfield(data, 'feature_groups_part')
        error('File does not contain "feature_groups_part": %s', ...
            filename);
    end

    % feature_groups_part = data.feature_groups_part;

  
        if isempty(feature_groups_all)
            feature_groups_all = data.feature_groups_part;
        else

        feature_groups_all = cat(1, feature_groups_all,data.feature_groups_part);
        end 

    clear data feature_groups_part
end

% -------------------------------------------------------------
% Finished
% -------------------------------------------------------------

fprintf('\nSuccessfully loaded all data.\n');
fprintf('Number of files: %d\n', numel(files));
fprintf('Number of samples: %d\n', numel(feature_groups_all));

end

