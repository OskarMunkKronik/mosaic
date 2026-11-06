function [volume_manual, apexTables] = extractTargetVolumes(Options, matFile_list, ...
        mzTarget, tol, rt_expected, rt_tol, compound_name, kList, opts)
%EXTRACTTARGETVOLUMES  Integrate target compounds in GCxGC sparse tensors,
%   plot them, and export apex / integration borders per sample to Excel.
%
%   [volume_manual, apexTables] = extractTargetVolumes(Options, matFile_list, ...
%        mzTarget, tol, rt_expected, rt_tol, compound_name, kList)
%
%   Inputs
%     Options        struct with Options.Paths.save2mat
%     matFile_list   dir() struct of the sample files
%     mzTarget       [nTargets x 1] target m/z
%     tol            m/z tolerance
%     rt_expected    [nTargets x 2] expected [1tR (min), 2tR (s)]
%     rt_tol         [1 x 2] half-window [1tR, 2tR] around rt_expected
%     compound_name  cellstr, nTargets
%     kList          sample indices to process (default: all)
%
%   Name-value options
%     MakePlots      (true)   save chromatogram .tiff per sample
%     SaveFig        (false)  also save .fig (slow, large)
%     VolumeManual   ([])     existing [nTargets x nSamples] matrix to update
%
%   volume_manual is returned as [nTargets x nSamples].
%
%   Excel output (one file per sample, in <save2mat>/target_compounds/rt_borders):
%     Compound | mzTarget | RT1_apex | RT2_apex | RT1_min | RT1_max |
%     RT2_min  | RT2_max  | Volume
%   RT1 in the units of Options.coords.X (min), RT2 in units of Options.coords.Y (s).
%   Edit the borders, then run reprocessVolumesFromExcel to update volume_manual.

arguments
    Options (1,1) struct
    matFile_list
    mzTarget (:,1) double
    tol (:,1) double
    rt_expected (:,2) double
    rt_tol (:,2) double
    compound_name
    kList (1,:) double = 1:numel(matFile_list)
    opts.MakePlots (1,1) logical = true
    opts.SaveFig (1,1) logical = false
    opts.VolumeManual double = []
end

nTargets = numel(mzTarget);
volume_manual = opts.VolumeManual;
if isempty(volume_manual)
    volume_manual = nan(nTargets, numel(matFile_list));
end
apexTables = cell(numel(matFile_list), 1);
allT = {};
root    = fullfile(Options.Paths.save2mat, 'target_compounds');
dirChrom = fullfile(root, 'chromatograms');
dirXls   = fullfile(root, 'rt_borders');
for p = {dirChrom, dirXls}
    if ~exist(p{1}, 'dir'), mkdir(p{1}); end
end

compound_name = string(compound_name(:));
nr = ceil(sqrt(nTargets));
nc = ceil(nTargets / nr);

for k = kList
    matPath = fullfile(Options.Paths.save2mat, 'SparseTensor', ...
        [matFile_list(k).name(1:end-4), '.mat']);
    % Only load what is needed (uRt, mod, fileList were never used)
    d = load(matPath, 'Z', 'mzroi_aug', 'Options');

    % Axis vectors (do this once per sample instead of indexing full grids)
    rt_1 = d.Options.coords.X(1, :);      % 1tR, row vector
    rt_2 = d.Options.coords.Y(:, 1).';    % 2tR, row vector
    mz   = d.mzroi_aug(:);

    apex   = nan(nTargets, 2);
    border = nan(nTargets, 4);            % [rt1min rt1max rt2min rt2max]

    if opts.MakePlots
        fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0 0 1 1]);
        tl  = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, matFile_list(k).name(1:end-4), 'Interpreter', 'none');
    end

    for j = 1:nTargets
        % find() gives integer index vectors -> cheaper, contiguous slicing
        i1 = find(abs(rt_1 - rt_expected(j,1)) <= rt_tol(j,1));
        i2 = find(abs(rt_2 - rt_expected(j,2)) <= rt_tol(j,2));
        im = find(abs(mz - mzTarget(j)) < tol(j));

        if opts.MakePlots
            ax = nexttile(tl);
        end
        if isempty(i1) || isempty(i2) || isempty(im)
            if opts.MakePlots, title(ax, compound_name(j) + " (no data)"); end
            continue
        end

        rt1w = rt_1(i1);
        rt2w = rt_2(i2);
        % Single slice; sum over m/z bins within tolerance
        Zs = sum(double(d.Z(i1, i2, im)), 3);

        volume_manual(j, k) = sum(Zs, 'all');
        border(j, :) = [min(rt1w) max(rt1w) min(rt2w) max(rt2w)];

        [pk, lin] = max(Zs(:));
        if pk > 0
            [r, c] = ind2sub(size(Zs), lin);
            apex(j, :) = [rt1w(r), rt2w(c)];
        end

        if opts.MakePlots && numel(rt1w) > 1 && numel(rt2w) > 1
            surf(ax, rt1w, rt2w, Zs.', 'FaceColor', 'interp', 'EdgeColor', 'none');
            view(ax, [0 90]); axis(ax, 'tight'); hold(ax, 'on');
            if pk > 0
                scatter3(ax, rt_expected(j,1), rt_expected(j,2), pk, 'red', 'filled');
            end
            xlabel(ax, '^1t_R (min)'); ylabel(ax, '^2t_R (s)');
            title(ax, compound_name(j), 'Interpreter', 'none');
        end
    end

    % ---- Excel export -------------------------------------------------
    T = table(compound_name, mzTarget, apex(:,1), apex(:,2), ...
        border(:,1), border(:,2), border(:,3), border(:,4), volume_manual(:, k), ...
        'VariableNames', {'Compound','mzTarget','RT1_apex','RT2_apex', ...
        'RT1_min','RT1_max','RT2_min','RT2_max','Volume'});
    apexTables{k} = T;
    allT{end+1} = addvars(T, repmat(string(matFile_list(k).name), nTargets, 1), ...
        'Before', 1, 'NewVariableNames', 'Sample'); %#ok<AGROW>
    xlsPath = fullfile(dirXls, [matFile_list(k).name, '.xlsx']);
    if exist(xlsPath, 'file'), delete(xlsPath); end
    writetable(T, xlsPath, 'Sheet', 'targets');

    % ---- Figures ------------------------------------------------------
    if opts.MakePlots
        exportgraphics(fig, fullfile(dirChrom, [matFile_list(k).name, '.tiff']), 'Resolution', 150);
        if opts.SaveFig
            savefig(fig, fullfile(dirChrom, [matFile_list(k).name, '.fig']));
        end
        close(fig);
    end
end

% Combined workbook for all samples in kList
writeSummary(allT, volume_manual, kList, {matFile_list.name}, ...
    fullfile(root, 'all_samples_rt_borders.xlsx'));
end

function writeSummary(allT, volume_manual, kList, sample_names, xlsPath)
% One workbook for the whole kList:
%   'targets' : long table (Sample + all per-compound columns)
%   'volumes' : compounds x samples matrix
if isempty(allT), return, end
Tall = vertcat(allT{:});
names = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(string(sample_names(kList))));
Tvol = [table(allT{1}.Compound, 'VariableNames', {'Compound'}), ...
        array2table(volume_manual(1:height(allT{1}), kList), 'VariableNames', cellstr(names))];
if exist(xlsPath, 'file'), delete(xlsPath); end
writetable(Tall, xlsPath, 'Sheet', 'targets');
writetable(Tvol, xlsPath, 'Sheet', 'volumes');
end
