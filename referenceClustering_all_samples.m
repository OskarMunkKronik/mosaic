
function clusters = referenceClustering_all_samples(points,Options)
% function clusters = referenceClustering_all_samples(points,Options)

n = length(points);
clusters = ones(n,1);

% --- Preallocate refs ---
refs = zeros(1,n);     % max possible size
ref_count = 1;

% Initialize first reference
if ~isfield(Options,'ref_sample')
    refs(1) = find(median([points.sample]) == [points.sample],1,'first');
else
    refs(1) = 1; % fallback (adjust if needed)
end

Rt_mat = [points.Rt];
max_dist = Options.Clustering.distance_componentization';
cos_thresh = 1-Options.Clustering.cutoff;

i_vec = 1:n; % simpler & safer

% --- Precompute ---
elution_profile = cell(2,1);
[elution_profile{1}, ~, ~] = align_elution_profiles_all([points.elution1]);
[elution_profile{2}, ~, ~] = align_elution_profiles_all([points.elution2]);

MS_norm = vecnorm([points.massSpec], 2, 1);
MS_all  = [points.massSpec];

for i = i_vec

    % Only use active refs
    current_refs = refs(1:ref_count);

    if ismember(i, current_refs)
        continue
    end

    % nTotal = n;

    if mod(i, max(1, round(0.05 * n))) == 0 || i == 1 || i == n
        fprintf(1, 'Grouping: %i/%i (%.0f%%)\n', i, n, 100*i/n);
    end

    assigned = false;

    % --- Distance filter ---
    validRefs = all(abs(Rt_mat(:,i) - Rt_mat(:,current_refs)) < max_dist, 1);
   
    if any(validRefs)

        c_vec = find(validRefs);
        ref_idx = current_refs(c_vec);

        % --- Elution cosine ---
        cos_elu = cell(2,1);

        for dim = 1:2
            X = elution_profile{dim};

            cand = X(:, i);
            refsX = X(:, ref_idx);

            cand_norm = norm(cand);
            refs_norm = vecnorm(refsX, 2, 1);

            cos_elu{dim} = (cand' * refsX) ./ (cand_norm * refs_norm + eps);
        end

        % --- MS cosine ---
        cand_MS = MS_all(:, i);
        ref_MS  = MS_all(:, ref_idx);

        cand_norm = MS_norm(i);
        ref_norm  = MS_norm(ref_idx);

        cos_MS = (cand_MS' * ref_MS) ./ (cand_norm * ref_norm + eps);

        % --- Combined score ---
        cos_prod_score = cos_elu{1} .* cos_elu{2} .* cos_MS;

        [bestScore, bestIdx] = max(cos_prod_score);

        if bestScore > cos_thresh
            clusters(i) = c_vec(bestIdx);
            assigned = true;
        end
    end

    % --- Add new reference (NO dynamic resizing) ---
    if ~assigned
        ref_count = ref_count + 1;
        refs(ref_count) = i;

        clusters(i) = ref_count;
    end

end

% Optional: trim refs (not required unless returned)
% refs = refs(1:ref_count);

end
%%last version 
