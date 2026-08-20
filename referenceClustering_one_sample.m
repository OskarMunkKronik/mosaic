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
    cos2 = elution2(:,i)' * elution2(:,refs(validRefs));
    % combined cosine similarity
    cosProd = cos1 .* cos2;%.*double(validRefs)
 
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
