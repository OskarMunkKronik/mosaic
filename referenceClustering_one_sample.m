function clusters = referenceClustering_one_sample(points,mz_ind,Z_feature,Options)
n = length(points);
clusters = ones(n,1);

refs = 1; % store reference indices
[~,idx] = sort([points.volume_summed],'descend');
refs = idx(1);
i_vec = idx(2:end);%2:n;
Rt_mat = [points.Rt];
max_dist = Options.Rtdev_within_sample';
cos_thresh = 1-Options.Clustering.cutoff;


% For all points, normalize once
% elu1 = cell2mat({points.elution1});
% elu2 = cell2mat({points.elution2});
elution1 = cellfun(@(x) x(:)/ (norm(x(:)) + eps), {points.elution1}, 'UniformOutput', false);
elution2 = cellfun(@(x) x(:)/ (norm(x(:)) + eps), {points.elution2}, 'UniformOutput', false);
elution1 = cell2mat(elution1);
elution2 = cell2mat(elution2);% size L x numRefs
for i = i_vec
    assigned = false;
    % Differences to all references at once

    % Keep only refs where both coordinates are within max_dist
    validRefs = all(abs(Rt_mat(:,i) - Rt_mat(:,refs)) <= max_dist, 1);
    if any(validRefs)
    pos_ValidRefs = find(validRefs);
    % Compute all dot products at once
    cos1 = elution1(:,i)' * elution1(:,refs(validRefs));     % 1 x numRefs %Fix so it only calculates the valid Refs!
    %     cos1 = elution1(:,i)' * elution1(:,refs);
    % Compute all dot products at once
    %     cos2 = elution2(:,i)' * elution2(:,refs);     % 1 x numRefs
    cos2 = elution2(:,i)' * elution2(:,refs(validRefs));
    % combined cosine similarity
    cosProd = cos1 .* cos2;%.*double(validRefs)
    % % mz_cand = mz_ind(i);
    % % mz_ref =  mz_ind([refs(validRefs)]);
    % % % [val,a] = min(abs(Options.mz_prob.mz_delta-abs(diff(mz_vec))));
    % % if sum(mz_ref==mz_cand) < length(mz_ref)
    % % mz_diff_mat = full(abs(mz_cand-mz_ref'));
    % % % n = 2;
    % % % nn = 1;
    % % 
    % % prob_mat = zeros(1,length(mz_ref)); % Initialize probability matrix
    % % ind  = find(mz_ref~=mz_cand);
    % % for n = ind;%1:numel(mz_ref)
    % % 
    % %     [val,a] = min(abs(Options.mz_prob.mz_delta - mz_diff_mat(n)));
    % % 
    % %     if val/mz_cand * 1e6 < Options.ppm_dev
    % %         prob_mat(n) =  Options.mz_prob.probability(a);
    % %     end
    % % end
    % % else 
    % %     prob_mat = ones(1,length(mz_ref));
    % % end 
    % prob_mat = compute_mz_prob(mz_ind(i),  mz_ind, refs, validRefs, Options);
    % [val,passIdx] = max(cosProd.*double(prob_mat>Options.mz_prob.theshold));
    % ind_max = find(clusters(refs(validRefs)) == clusters);
    % [~,ref_max] = max(points(ind_max).volume_summed);
    % clusters

    % [ratio] = compute_ratio_cosine(elution1,elution2, i,refs(validRefs));
   
    % [ratio] = compute_ratio(elu1,elu2, i,refs(validRefs),clusters);
    [val,passIdx] = max(cosProd);%.*ratio);
    %%
    passIdx = pos_ValidRefs(passIdx);
    % find refs that pass threshold
    if val>cos_thresh%~isempty(passIdx)
        clusters(i) = passIdx; % or index in refs depending on logic
        assigned = true;
    end
    end 
    if ~assigned
        refs(end+1) = i; % new cluster reference
        clusters(i) = length(refs);
        i_vec = i_vec(~(ismember(i_vec,1:length(refs))));
    end
end

end
%%
% function prob_mat = compute_mz_prob(mz_cand, mz_ind, refs, validRefs, Options)
% %COMPUTE_MZ_PROB Compute probability matrix for m/z candidate vs. references.
% %
% %   prob_mat = compute_mz_prob(mz_cand, mz_ind, refs, validRefs, Options)
% %
% %   Inputs:
% %       mz_cand   - scalar m/z candidate value
% %       mz_ind    - vector of all m/z indices
% %       refs      - reference index vector
% %       validRefs - logical or numeric vector of valid reference indices
% %       Options   - struct with fields:
% %                       .mz_prob.mz_delta      : array of delta m/z values
% %                       .mz_prob.probability   : corresponding probabilities
% %                       .ppm_dev               : maximum ppm deviation allowed
% %
% %   Output:
% %       prob_mat  - row vector of probabilities (1 x Nref)
% %
% %   Notes:
% %       - Fully vectorized (no for-loops)
% %       - Assumes Options.mz_prob.mz_delta and Options.mz_prob.probability
% %         are of the same length
% %
% %   Example:
% %       prob = compute_mz_prob(500.123, mz_ind, refs, validRefs, Options);
% 
%     % Extract reference m/z values
%     mz_ref = mz_ind(refs(validRefs));
% 
%     % If all references equal the candidate, return ones
%     if all(mz_ref == mz_cand)
%         prob_mat = ones(1, numel(mz_ref));
%         return;
%     end
% 
%     % Compute absolute m/z differences
%     mz_diff = abs(mz_cand - mz_ref(:));  % column vector
% 
%     % Find closest mz_delta for each difference
%     % [val, idx] = min(abs(Options.mz_prob.mz_delta(:)' - mz_diff), [], 2);
%     [r,c] = find((abs(Options.mz_prob.mz_delta(:)' - mz_diff)<(1+Options.ppm_dev*1e-6)*mz_diff) & (abs(Options.mz_prob.mz_delta(:)' - mz_diff)>(1-Options.ppm_dev*1e-6)*mz_diff));
%     prob_mat =accumarray(r , Options.mz_prob.probability(c),[],@max);
%     % max(Options.mz_prob.probability(c(r==1)))
%     % Calculate ppm deviation
%     % ppm_dev_vec = val ./ mz_cand * 1e6;
% 
%     % Initialize output
%     % prob_mat = zeros(1, numel(mz_ref));
% 
%     % Logical index of valid ppm matches
%     % valid_idx = ppm_dev_vec < Options.ppm_dev;
% 
%     % Assign corresponding probabilities
%     prob_mat(valid_idx) = Options.mz_prob.probability(idx(valid_idx));
% 
% end
function prob_mat = compute_mz_prob(mz_cand, mz_ind, refs, validRefs, Options)
%COMPUTE_MZ_PROB Compute probability matrix for m/z candidate vs. references.
%
%   prob_mat = compute_mz_prob(mz_cand, mz_ind, refs, validRefs, Options)
%
%   Inputs:
%       mz_cand   - scalar m/z candidate value
%       mz_ind    - vector of all m/z indices
%       refs      - reference index vector
%       validRefs - logical or numeric vector of valid reference indices
%       Options   - struct with fields:
%                       .mz_prob.mz_delta      : array of delta m/z values
%                       .mz_prob.probability   : corresponding probabilities
%                       .ppm_dev               : maximum ppm deviation (ppm)
%
%   Output:
%       prob_mat  - row vector of probabilities (1 x Nref)
%
%   Notes:
%       - Fully vectorized (no loops)
%       - Takes the maximum probability among all deltas within ppm tolerance

    % Extract valid reference m/z values
    mz_ref = mz_ind(refs(validRefs));
    nRef = numel(mz_ref);

    % If all references are identical to the candidate, return ones
    if all(mz_ref == mz_cand)
        prob_mat = ones(1, nRef);
        return;
    end

    % Compute absolute m/z differences
    mz_diff = abs(mz_cand - mz_ref(:));  % nRef x 1

    % Expand delta grid
    mz_delta = Options.mz_prob.mz_delta(:)';  % 1 x Ndelta
    tol = Options.ppm_dev * 1e-6;             % fractional tolerance

    % Logical matrix of matches within tolerance
    within_tol = abs(mz_delta - mz_diff)./min(mz_cand,mz_ref(:)) <= tol * mz_diff;

    % If nothing matches, return zeros
    if ~any(within_tol, 'all')
        prob_mat = zeros(1, nRef);
        return;
    end

    % Compute per-row (reference) max probability within tolerance
    prob_table = Options.mz_prob.probability(:)';  % ensure row
    prob_mat_col = max(within_tol .* prob_table, [], 2);  % nRef x 1

    % Convert NaNs (no matches) to zero
    prob_mat_col(isnan(prob_mat_col)) = 0;

    % Ensure output is a row vector
    prob_mat = prob_mat_col.';
end
