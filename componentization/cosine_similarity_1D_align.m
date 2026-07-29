function [x_ref, x_cand, cand_aligned, cosine_sim] = cosine_similarity_1D_align(ref, cand)
%COSINE_SIMILARITY_1D Compute cosine similarity between two 1D elution profiles.
%   Aligns candidate (cand) to reference (ref) by peak apex before comparison.

    % --- Normalize inputs ---
    ref = ref(:);
    cand = cand(:);
    ref = ref / (max(ref) + eps);
    cand = cand / (max(cand) + eps);

    % --- Find apex positions ---
    [~, idx_ref] = max(ref);
    [~, idx_cand] = max(cand);

    % --- Compute alignment shift (index difference) ---
    shift = idx_ref - idx_cand;

    % --- Align candidate signal safely ---
    if shift > 0
        cand_aligned = [zeros(shift,1); cand]; % pad front
    elseif shift < 0
        cand_aligned = cand(-shift+1:end); % trim front
    else
        cand_aligned = cand;
    end

    % --- Adjust to same length as ref ---
    if length(cand_aligned) < length(ref)
        cand_aligned(end+1:length(ref)) = 0; % pad end
    elseif length(cand_aligned) > length(ref)
        cand_aligned = cand_aligned(1:length(ref)); % trim extra
    end

    % --- Define normalized x-axes ---
    x_ref  = linspace(0, 1, length(ref));
    x_cand = linspace(0, 1, length(cand_aligned));

    % --- Interpolate candidate if slightly different ---
    if length(cand_aligned) ~= length(ref)
        cand_aligned = interp1(x_cand, cand_aligned, x_ref, 'linear', 'extrap');
    end

    % --- Compute cosine similarity ---
    cosine_sim = dot(cand_aligned, ref) / (norm(cand_aligned) * norm(ref) + eps);

end
