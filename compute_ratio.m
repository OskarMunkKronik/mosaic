function [ratio] =  compute_ratio(elution1,elution2, i,refs,clusters)

% ref_idx = vec(1);
% nVec = numel(vec);
nDim = 2;

% Extract all reference and candidate signals into matrices

% ref_data = cellfun(@(d1, d2) [d1(:), d2(:)], ...
%     sum_elution_profile{k}{ref_idx}(1), sum_elution_profile{k}{ref_idx}(2), 'UniformOutput', false);
% ref_data = [sum_elution_profile{k}{ref_idx}{1}(:), sum_elution_profile{k}{ref_idx}{2}(:)];
[ref_in_cluster,~] = find(clusters(refs)' == clusters);
[r,c,val] = find(elution1(:,ref_in_cluster));
ref_data{1} = accumarray([r,c],val,[size(elution1,1),max(c)],@max);

% ref_in_cluster = find(clusters(refs) == clusters);
[r,c,val] = find(elution2(:,ref_in_cluster));
ref_data{2} = accumarray([r,c],val,[size(elution2,1),max(c)],@max);

% ref_data{1} = accumarray(repmfind(clusters(refs) == clusters),elution1(:,clusters(refs) == clusters));
% ref_data{2} = elution2(:,refs);

cand_data{1} = elution1(:,i);
cand_data{2} = elution2(:,i);
for dim = 1:nDim
ref_data{dim}(ref_data{dim}==0)=NaN;
cand_data{dim}(cand_data{dim}==0)=NaN;
end 

% % cand_data = cellfun(@(v) ...
% %     [sum_elution_profile{k}{v}{1}(:), sum_elution_profile{k}{v}{2}(:)], ...
% %     num2cell(vec), 'UniformOutput', false);
% 
% % Pad all candidates to equal length
% maxLen = max(cellfun(@(c) size(c,1), cand_data));
% cand_mat = NaN(maxLen, nDim, nVec);
% for i = 1:nVec
%     ci = cand_data{i};
%     cand_mat(1:size(ci,1), :, i) = ci;
% end
% cand_mat(cand_mat==0)=NaN;

% % Replicate reference for all candidates
% ref_mat = repmat(ref_data, [1 1 nVec]);

% Compute ratio > 1 in a fully vectorized manner
ratio = nan(nDim,length(refs));
for n = 1:length(refs)
    for dim = 1:nDim
ind = ~isnan(ref_data{dim}(:,n)) & ~isnan(cand_data{dim});
ratio_mask = ref_data{dim}(:,n) ./ cand_data{dim} > 1;
ratio(dim,n) = length(unique(ratio_mask(ind)))';
    end
end 
% clf
% for n = 1:length(refs)
%     for dim = 1:nDim
%          ind = ~isnan(ref_data{dim}(:,n)) & ~isnan(cand_data{dim});
%         subplot(2,1,dim)
%         plot(ref_data{dim}(ind,n))
%         hold on 
%         plot(cand_data{dim}(ind),'b')
%         yyaxis right
% 
%         rat = ref_data{dim}(:,n) ./ cand_data{dim};
%         plot(rat(ind))
%         yline(1)
%     end 
% end 


ratio = double(sum(ratio,1)  == 2);
% % Compute cosine similarity fully vectorized
% ref_norm = sqrt(nansum(ref_data.^2, 1));
% cand_norms = sqrt(nansum(cand_mat.^2, 1, 'omitnan'));
% dot_prods = squeeze(nansum(cand_mat .* ref_mat, 1, 'omitnan'));
% 
% cos_sim = dot_prods ./ (ref_norm .* cand_norms + eps);
% cos_sim = squeeze(cos_sim)';
end

