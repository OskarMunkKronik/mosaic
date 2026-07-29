function [W, H, cg, exp_var, X_unfolded] = perform_nnmf_unfold(Z_tmp, Options)
%PERFORM_NNMF_UNFOLD Unfold 3D elution data and run iterative NNMF
%
%   [W, H, cg, exp_var, X_unfolded] = perform_nnmf_unfold(Z_tmp, Options)
%
%   Inputs:
%       Z_tmp   - 3D array (X × Y × N) of aligned elution maps
%       Options - struct with fields:
%                   .doPlot : (true/false) plot unfolded signals
%                   .curve_resolution.exp_var_thresh : target explained variance (0–1)
%                   .curve_resolution.algorithm : NNMF algorithm ('mult', 'als', etc.)
%
%   Outputs:
%       W, H        - NNMF factorization matrices
%       cg          - index of max component for each feature
%       exp_var     - final explained variance
%       X_unfolded  - unfolded (2D) data matrix [XY × N]

% -----------------------------
% Step 1: Unfold 3D cube into 2D
% -----------------------------
if length(size(Z_tmp))>3
    X_unfolded = zeros(size(Z_tmp,1) * size(Z_tmp,2)* size(Z_tmp,4), size(Z_tmp,3));
    for n = 1:size(Z_tmp,3)
        X_tmp = permute(squeeze(Z_tmp(:,:,n,:)),[2,1,3]);
        X_unfolded(:,n) = X_tmp(:);
    end
else
    X_unfolded = zeros(size(Z_tmp,1) * size(Z_tmp,2), size(Z_tmp,3));
    for n = 1:size(Z_tmp,3)
        X_tmp = Z_tmp(:,:,n)';
        X_unfolded(:,n) = X_tmp(:);
    end
end



% X_unfolded = X_unfolded./max(X_unfolded,[],1);
% -----------------------------
% Step 2: Optional plotting
% -----------------------------
if isfield(Options, 'doPlot') && Options.doPlot
    clf
    hold on
    plot(X_unfolded)
    hold off
end

% -----------------------------
% Step 3: Iterative NNMF
% -----------------------------
exp_var = 0;
F = 0;
max_factor = min([Options.curve_resolution.max_factor ,size(X_unfolded)]);
while exp_var < Options.curve_resolution.exp_var_thresh && F < max_factor%Options.curve_resolution.max_factor 
    F = F + 1;
    warning('off','all')
    [W, H] = nnmf(X_unfolded, F, 'algorithm', Options.curve_resolution.algorithm);

    SSQ_error = sum((X_unfolded - W*H).^2, 'all');
    SSQ_total = sum(X_unfolded.^2, 'all');
    exp_var = 1 - (SSQ_error / SSQ_total);
end

% -----------------------------
% Step 4: Assign dominant component
% -----------------------------
[~, cg] = max(H, [], 1);

end
