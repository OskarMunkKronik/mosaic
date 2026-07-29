function clusters = referenceClustering_all_samples_new(points,Options)
n = length(points);
clusters = ones(n,1);
% clusters(1,1) = 1;
% refs = 1; % store reference indices
if ~isfield(Options,'ref_sample')
    refs = find(median([points.sample]) == [points.sample],1,'first');
end
% clusters = 1;
Rt_mat = [points.Rt];
max_dist = Options.Clustering.distance_componentization';
cos_thresh = 1-Options.Clustering.cutoff;
doPlot = Options.doPlot;
i_vec = find(refs ~= 1:size(Rt_mat,2));

% align elution profiles
elution_profile = cell(2,1);
[elution_profile{1}, ~, ~] = align_elution_profiles_all([points.elution1]);
[elution_profile{2}, ~, ~] = align_elution_profiles_all([points.elution2]);
MS_norm = vecnorm([points.massSpec], 2, 1);

for i = i_vec
    assigned = false;
    % Keep only refs where both coordinates are within max_dist
    validRefs = all(abs(Rt_mat(:,i) - Rt_mat(:,refs)) < max_dist, 1);
    c_vec = find(validRefs);
    %%%%
    %%%%
    % cosineSim(elu_profile{1},   points(refs(c)).elution1)
    % cosineSim(elution_profile{1}(:,i),  elution_profile{1}(:,refs(c_vec)))
    cos_elu = cell(2,1);
    if ~isempty(c_vec)
        for dim = 1:2
            X = elution_profile{dim};      % n × D

            cand = X(:, i);               % n × 1
            refsX = X(:, refs(c_vec));    % n × K

            cand_norm = norm(cand);
            refs_norm = vecnorm(refsX, 2, 1);

            cos_elu{dim} = (cand' * refsX) ./ (cand_norm * refs_norm + eps);
            % cosineSim(elution_profile{1}(:,i),  elution_profile{1}(:,refs(c_vec)))
            %     cos_elu{dim} = arrayfun( ...
            %     @(c) cosineSim(elution_profile{dim}(:, i), ...
            %                    elution_profile{dim}(:, refs(c))), ...
            %     c_vec ...
            % );
            % isequal(cs,cos_elu{dim})
        end
        % Stack once
        % Keep everything sparse
        ref_MS  = [points(refs(c_vec)).massSpec];   % n × K sparse
        cand_MS = points(i).massSpec;               % n × 1 sparse

        % Precompute norms (sparse-safe)
        cand_norm = MS_norm(i);
        ref_norm  =MS_norm(refs(c_vec));

        % Cosine similarity
        cos_MS = (cand_MS' * ref_MS) ./ (cand_norm * ref_norm + eps);
        % cos_MS = zeros(length(c_vec),1);
        % for k = 1:numel(c_vec)
        %     cos_MS(k) = cosineSim(points(i).massSpec, ...
        %                           points(refs(c_vec(k))).massSpec);
        % end


        cos_prod_score(c_vec)= cos_elu{1}.*cos_elu{2}.*cos_MS;
    else
        cos_prod_score = [];

    end
    for c_count = 1:length(c_vec)%length(refs)
        % profile_ref{1} = points(refs(c)).elution1;
        %            profile_cmp{1} = points(i).elution1;
        %            profile_ref{2} = points(refs(c)).elution2;
        %            profile_cmp{2} = points(i).elution2;
        %
        %            [elu_profile, ~] = align_elution_profiles(profile_ref, profile_cmp, doPlot,Options);

        c = c_vec(c_count);
        if validRefs(c)
            % profile_ref{1} = points(refs(c)).elution1;
            % profile_cmp{1} = points(i).elution1;
            % profile_ref{2} = points(refs(c)).elution2;
            % profile_cmp{2} = points(i).elution2;
            %
            % [elu_profile, ~] = align_elution_profiles(profile_ref, profile_cmp, doPlot,Options);
            if cos_prod_score(c)  >cos_thresh
                clusters(i) = c;
                assigned = true;
                break;
            end
        end
    end
    if ~assigned
        refs(end+1) = i; % new cluster reference
        clusters(i) = length(refs);
        i_vec = i_vec(~(ismember(i_vec,1:length(refs))));
    end
end



end

function c = cosineSim(v1,v2)
c = dot(v1,v2)/(norm(v1)*norm(v2)+eps);
end
% function clusters = referenceClustering_all_samples_new(points, Options)
% 
% n = numel(points);
% clusters = ones(n,1);
% 
% % Initial reference
% if ~isfield(Options,'ref_sample')
%     refs = find(median([points.sample]) == [points.sample], 1, 'first');
% end
% 
% % Precompute constants
% Rt_mat    = [points.Rt];
% max_dist  = Options.Clustering.distance_componentization';
% cos_thresh = 1 - Options.Clustering.cutoff;
% 
% % Indices to process
% i_vec = find(refs ~= 1:size(Rt_mat,2));
% 
% % Pre-align elution profiles (global alignment)
% elution_profile = cell(2,1);
% [elution_profile{1}, ~, ~] = align_elution_profiles_all([points.elution1]);
% [elution_profile{2}, ~, ~] = align_elution_profiles_all([points.elution2]);
% 
% % Precompute MS norms (sparse-safe)
% MS_norm = vecnorm([points.massSpec], 2, 1);
% 
% % ==============================
% % Main clustering loop
% % ==============================
% for i = i_vec
% 
%     assigned = false;
% 
%     % Candidate references within RT tolerance
%     validRefs = all(abs(Rt_mat(:,i) - Rt_mat(:,refs)) < max_dist, 1);
%     c_vec = find(validRefs);
% 
%     if ~isempty(c_vec)
% 
%         % ---- Elution cosine similarity (2 dimensions) ----
%         cos_elu = cell(2,1);
% 
%         for dim = 1:2
%             X = elution_profile{dim};
% 
%             cand  = X(:, i);
%             refsX = X(:, refs(c_vec));
% 
%             cand_norm = norm(cand);
%             refs_norm = vecnorm(refsX, 2, 1);
% 
%             cos_elu{dim} = (cand' * refsX) ./ ...
%                            (cand_norm * refs_norm + eps);
%         end
% 
%         % ---- Mass spectrum cosine similarity (sparse) ----
%         cand_MS = points(i).massSpec;
%         ref_MS  = [points(refs(c_vec)).massSpec];
% 
%         cos_MS = (cand_MS' * ref_MS) ./ ...
%                  (MS_norm(i) * MS_norm(refs(c_vec)) + eps);
% 
%         % ---- Combined similarity score ----
%         cos_prod_score = cos_elu{1} .* cos_elu{2} .* cos_MS';
% 
%     else
%         cos_prod_score = 0;
%     end
% 
%     % ---- Assign to first matching reference ----
%     for k = 1:numel(c_vec)
%         if cos_prod_score(k) > cos_thresh
%             clusters(i) = c_vec(k);
%             assigned = true;
%             break;
%         end
%     end
% 
%     % ---- Create new reference if unassigned ----
%     if ~assigned
%         refs(end+1) = i;
%         clusters(i) = numel(refs);
% 
%         % Remove newly assigned ref from future candidates
%         i_vec = i_vec(~ismember(i_vec, 1:numel(refs)));
%     end
% end
% 
% end
