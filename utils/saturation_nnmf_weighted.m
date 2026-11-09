%% Weighted NMF with saturation mask on real GCxGC-HRMS data
% Base MATLAB only. Run as a script (local functions are at the end).
clear; close all; rng(1);

%% 0. Settings
dataDir  = "N:\SCIENCE-PLEN-ECP_AnalytChem\Projects\MST\test";
tmplFile = "D:\MST\integration\template\target_integration_template_BatchN.xlsx";
n = 5;            % index of the .mat file
j = 45;           % index of the target compound in the template
k = 2;            % number of NMF components (check with svd/scree, see notes)
nIter   = 1000;   % multiplicative-update iterations
nStarts = 5;      % random restarts, best weighted cost is kept

% Mask: use only points between loFrac and upFrac of the reference max
loFrac  = 0.10;
upFrac  = 0.90;
refMode = 'fixed';   % 'perIon' : max of each EIC in the window
                      % 'run'    : max of the whole file (d.Z)
                      % 'fixed'  : tClip below
tClip   = [];         % detector clip level, only used for refMode = 'fixed'

%% 1. Read template and select target
fileList   = dir(fullfile(dataDir, "*.mat"));
input_data = readtable(tmplFile);

rt_expected  = input_data{:, ["x1stDimensionTime_min_","x2ndDimensionTime_s_"]};
rtTol        = input_data{:, ["x1stDimensionR_T_Deviation_s_","x2ndDimensionR_T_Deviation_s_"]} ./ [60, 1];  % [1tR (min), 2tR (s)]
mzTarget     = input_data{:, "mzTarget"};
mzTol        = input_data{:, "mzTol"};
compoundName = input_data.Name;

%% 2. Load file and cut out the RT window
d = load(fullfile(fileList(n).folder, fileList(n).name));

rt_1 = d.Options.coords.X(1, :);
rt_2 = d.Options.coords.Y(:, 1);
ind_1  = find(abs(rt_1 - rt_expected(j,1)) <= rtTol(j,1));
ind_2  = find(abs(rt_2 - rt_expected(j,2)) <= rtTol(j,2));
mz_ind = find(abs(d.mzroi_aug - mzTarget(j)) < mzTol(j));
assert(~isempty(ind_1) && ~isempty(ind_2), 'No scans inside the RT window.');
assert(~isempty(mz_ind), 'Target m/z not found in d.mzroi_aug.');

% Orientation check: coords.X varies along columns, coords.Y along rows,
% so Z is expected as [2tR x 1tR x mz (x extra)]. Handle the other layout too.
fprintf('size(Z) = %s, numel(rt_2) = %d, numel(rt_1) = %d\n', ...
        mat2str(size(d.Z)), numel(rt_2), numel(rt_1));
Z_tmp = double(d.Z(ind_1, ind_2, :));
if size(d.Z,1) == numel(rt_2) && size(d.Z,2) == numel(rt_1)
    Z_tmp = double(Z_tmp);
else
    Z_tmp = permute(Z_tmp, [2 1 3 ]);
end

% Unfold to scans x ions, in chronological order (2tR fastest, then 1tR, then 4th dim)
nMz   = size(Z_tmp, 3);
Xsat  = reshape(permute(Z_tmp, [1 2 4 3]), [], nMz);
Xsat  = max(Xsat, 0);
mzAxis = d.mzroi_aug(:)';

nScans = size(Xsat, 1);  nIons = size(Xsat, 2);
t = (1:nScans)';
fprintf('Window: %d scans x %d ions\n', nScans, nIons);

%% 3. Window mask (10 % - 90 % of the reference max)
tClip = 1.5e7%max(Z_tmp,[],'all')*0.5
switch refMode
    case 'perIon'
        refMax = max(Xsat, [], 1);
    case 'run'
        refMax = max(double(d.Z), [], 'all') * ones(1, nIons);
    case 'fixed'
        assert(~isempty(tClip), 'Set tClip for refMode = ''fixed''.');
        refMax = tClip * ones(1, nIons);
    otherwise
        error('Unknown refMode.');
end
tLo = 1 * min(Z_tmp,[],'all');
tUp = upFrac * refMax;

tooHigh = Xsat >= tUp;
tooHigh = tooHigh | circshift(tooHigh, 1, 1) | circshift(tooHigh, -1, 1);  % dilate +/-1 scan
tooLow  = Xsat <= tLo;
bad = tooHigh | tooLow;
W   = double(~bad);                       % 1 = use, 0 = ignore
fprintf('Masked points: %d (%.1f%%)  [high: %d, low: %d]\n', ...
        sum(bad(:)), 100*mean(bad(:)), sum(tooHigh(:)), sum(tooLow(:)));

% 4. NMF: unweighted vs weighted, best of nStarts
best_u = struct('cost', Inf);  best_w = struct('cost', Inf);
for s = 1:nStarts
    A0 = rand(nScans, k) * mean(Xsat(:));
    B0 = rand(k, nIons);
    [Au, Bu, cu] = wnmf(Xsat, ones(size(Xsat)), A0, B0, nIter);
    [Aw, Bw, cw] = wnmf(Xsat, W,                A0, B0, nIter);
    if cu(end) < best_u.cost, best_u = struct('A',Au,'B',Bu,'cost',cu(end)); end
    if cw(end) < best_w.cost, best_w = struct('A',Aw,'B',Bw,'cost',cw(end)); end
end

% Normalise spectra to unit norm, order components by profile area
[A_u, B_u] = normalise(best_u.A, best_u.B);
[A_w, B_w] = normalise(best_w.A, best_w.B);

Xhat_u = A_u * B_u;
Xhat_w = A_w * B_w;
[~,mz_ind] = max(sum(Xsat))
%% 5. Target compound: EIC and area
eic_data = sum(Xsat(:, mz_ind), 2);
eic_u    = sum(Xhat_u(:, mz_ind), 2);
eic_w    = sum(Xhat_w(:, mz_ind), 2);
maskedT  = any(bad(:, mz_ind), 2);

fprintf('\nTarget: %s  (m/z %.4f, %d ion channel(s))\n', string(compoundName(j)), mzTarget(j), numel(mz_ind));
fprintf('  Area, raw data     : %.4g\n', sum(eic_data));
fprintf('  Area, NMF          : %.4g  (x%.2f of raw)\n', sum(eic_u), sum(eic_u)/sum(eic_data));
fprintf('  Area, weighted NMF : %.4g  (x%.2f of raw)\n', sum(eic_w), sum(eic_w)/sum(eic_data));

% 6. Plots
figure(1); clf; set(gcf, 'Position', [100 100 1200 750]);

subplot(2,2,1)
plot(t, eic_data, 'r.-'); hold on
plot(t(maskedT), eic_data(maskedT), 'ko', 'MarkerSize', 6);
plot(t, eic_u, 'b--', 'LineWidth', 1.5);
plot(t, eic_w, 'g-',  'LineWidth', 1.5);
yline(refMax)
legend('Data', 'Masked points', 'NMF', 'Weighted NMF', 'Location', 'best');
title(sprintf('Target EIC: %s', string(compoundName(j))), 'Interpreter', 'none');
xlabel('Scan (unfolded)'); ylabel('Intensity');

subplot(2,2,2)
imagesc(t, 1:nIons, ~bad'); colormap(gca, gray);
title('Weight mask (white = used)'); xlabel('Scan'); ylabel('Ion index');

subplot(2,2,3)
plot(t, A_w, 'LineWidth', 1.3);
title('Weighted NMF: recovered profiles'); xlabel('Scan'); ylabel('Intensity');
legend(compose('Comp %d', 1:k));

subplot(2,2,4)
bar(mzAxis, B_w(1,:)', 'k'); hold on
bar(mzAxis, -B_u(1,:)', 'b');
title('Component 1 spectrum: weighted (up) vs unweighted (down)');
xlabel('m/z'); ylabel('Unit-norm intensity');

%% ---------------- local functions ----------------
function [A, B, cost] = wnmf(X, W, A, B, nIter)
% Weighted NMF, min sum(W .* (X - A*B).^2), multiplicative updates
    epsv = 1e-12;
    WX = W .* X;
    cost = zeros(nIter, 1);
    for it = 1:nIter
        A = A .* (WX * B') ./ ((W .* (A*B)) * B' + epsv);
        B = B .* (A' * WX) ./ (A' * (W .* (A*B)) + epsv);
        cost(it) = sum(sum(W .* (X - A*B).^2));
    end
end

function [A, B] = normalise(A, B)
% Unit-norm spectra (rows of B), profiles rescaled so A*B is unchanged,
% components sorted by profile area (largest first)
    nb = vecnorm(B, 2, 2);          % k x 1
    B  = B ./ nb;
    A  = A .* nb';
    [~, o] = sort(sum(A, 1), 'descend');
    A = A(:, o);  B = B(o, :);
end