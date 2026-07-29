function [cosMat,cosMat_combined] = calculate_cos_matrix(sum_elution_profile,Z_feature)
n_feature = numel(sum_elution_profile);
% Preallocate similarity matrices
% cosMat = cell(1,2);
cosMat= zeros([n_feature,n_feature,2]); % square similarity matrix
cosMat_combined = zeros(n_feature);
%maybe do something like this  [r,c] = find(triu(ones(n_feature),1))
% [r,c] = find(triu(ones(n_feature),1));
% Compute cosine similarities
for nn = 1:n_feature
    for n = 1:n_feature
        for dim = 1:2
            x = sum_elution_profile{nn}{dim};
            y = sum_elution_profile{n}{dim};
            x =x./sum(x,'all');
            y = y./sum(y,'all');
            cosMat(n,nn,dim) = dot(x,y) / (norm(x)*norm(y) + eps);
        end
        cosMat_combined(n,nn) = cosMat(n,nn,1) * cosMat(n,nn,2);
    end
end
n_feature = numel(sum_elution_profile);


%     Precompute flattened and normalized vectors for each dim
X = cell(1,2);
for dim = 1:2
    %         Flatten each profile into a column vector
    M = cellfun(@(c) c{dim}(:), sum_elution_profile, 'UniformOutput', false);
    M = cat(2, M{:});  % columns = features

    %         Normalize by sum of elements (independently)
    M = M ./ sum(M,1);

    %         L2-normalize columns for cosine similarity
    norms = sqrt(sum(M.^2,1)) + eps;
    X{dim} = bsxfun(@rdivide, M, norms);
end

%     Compute cosine similarities for each dim
cosMat(:,:,1) = X{1}' * X{1};
cosMat(:,:,2) = X{2}' * X{2};

cosMat_combined = cosMat(:,:,1) .* cosMat(:,:,2);
end

%%
