%% Aux functions
function convert2Ztensor_and_plot(k, q, loadDir, fileList_mat, fileList, mod_k, uRt_k, mzroi_aug, dims, Options)
%PROCESSSAMPLE  BPC/TIC plots + sparse tensor for one sample (parfor-safe).
%
%   k             sample index
%   q             parallel.pool.DataQueue for progress (or [] to skip)
%   loadDir       folder with the MSroi_aug .mat files
%   fileList_mat  dir() struct of the .mat files
%   fileList      dir() struct of the raw files (used for output names)
%   mod_k         mod{k}   (modulation index per RT column)
%   uRt_k         uRt{k}   (1tR index per RT column)
%   mzroi_aug     m/z axis
%   dims          [nMod nRt], computed ONCE outside the loop:
%                 [max(cellfun(@max,mod(:))), max(cellfun(@max,uRt(:)))]
%   Options       struct with Options.firstRt, Options.coords, Options.Paths.save2mat

name = fileList(k).name(1:end-4);
root = Options.Paths.save2mat;

% ---- load data (only the variable needed) -------------------------------
d = load(fullfile(loadDir, fileList_mat(k).name), 'MSroi_aug');

% ---- map column indices to mod / RT coordinates -------------------------
[ionSlice, c, val] = find(d.MSroi_aug(:, Options.firstRt(k):end));
d = []; %#ok<NASGU>   free memory early

modCoord = mod_k(c);
rtCoord  = uRt_k(c);
sub      = [modCoord(:), rtCoord(:)];

% ---- BPC and TIC maps (transposed to match meshgrid X/Y) ----------------
BPC = accumarray(sub, val, dims, @max).';
TIC = accumarray(sub, val, dims).';          % default is sum, faster than @sum

plotMap(Options.coords.X, Options.coords.Y, BPC, ['BPC: ', name], ...
    fullfile(root, 'BPC', [name, '.tif']));
plotMap(Options.coords.X, Options.coords.Y, TIC, ['TIC: ', name], ...
    fullfile(root, 'TIC', [name, '.tif']));

% ---- sparse 3D tensor (requires Tensor Toolbox on the workers) ----------
Z = sptensor([modCoord(:), rtCoord(:), ionSlice(:)], val, [dims, numel(mzroi_aug)]); %#ok<NASGU>

% ---- save ---------------------------------------------------------------
save(fullfile(root, 'SparseTensor', [fileList(k).name, '.mat']), ...
    'Z', 'mzroi_aug', 'fileList', 'Options', '-v7.3');

if ~isempty(q), send(q, k); end
end

function plotMap(X, Y, M, ttl, outPath)
fig = figure('Visible', 'off');
ax  = axes(fig);
surf(ax, X, Y, M, 'EdgeColor', 'interp', 'FaceColor', 'interp');
view(ax, [0 90]); axis(ax, 'tight'); colorbar(ax);
xlabel(ax, '^1t_R (min)'); ylabel(ax, '^2t_R (s)');
title(ax, ttl, 'Interpreter', 'none');
exportgraphics(fig, outPath, 'Resolution', 150);
close(fig);
end
