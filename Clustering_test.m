
%%

i = 1;
points = struct();
count = 0;
for k_count = 1:length(sample_vec)
    k = sample_vec(k_count);

    for n = 1:length(clusters{k})

        if allPeaks{k}(n,6)>0
            points(i).Rt                      = allPeaks{k}(n,[2,3])';
            points(i).elution1                = sum_elution_profile{k}{n}{1};
            points(i).elution2                = sum_elution_profile{k}{n}{2};
            points(i).feature_elu_2D          = Z_feature{k}(n);
            points(i).massSpec                = mass_spectrum{k}(:,clusters{k}(n));
            points(i).volume_summed           = sum(allPeaks{k}(clusters{k}(n)==clusters{k},5));
            [~,a]                             = max(mass_spectrum{k}(:,clusters{k}(n)));
            points(i).basePeak                = mzroi_aug(a);
            points(i).sample                  = k;
            points(i).sample_feature_ID       = allPeaks{k}(n,6);
            points(i).sample_cluster_ID       = clusters{k}(n);
            
            if size(allPeaks{k},2)>6
            points(i).IS_flag                 = allPeaks{k}(n,7);
            end 
            
            i = i +1;
            
        end
    end
end
points
%%

% %--- Run clustering using reference sample ---
% clusters_test = referenceClustering(points,Options);
clusters_test = referenceClustering_all_samples(points,Options)
numClusters = max(clusters_test);
disp(['Number of clusters found: ', num2str(numClusters)]);

%%
% 


mzTarget = 275.2350%268.191
d = full([points.basePeak;clusters_test';points.Rt;points.sample]')
ind = find(abs(d(:,1)-mzTarget)<0.001 )
d(ind,:)
clf
col = rand()
col = rand(sum(ind),3);
scatter(d(ind,3),d(ind,4),[],d(ind,2))
text(d(ind,3),d(ind,4),num2str(d(ind,2)))

%%
cG_vec = [212]
clf
for S = 1:length(cG_vec)
    subplot(length(cG_vec),1,S)
a = find(d(:,2) == cG_vec(S));
for n = 1:length(a);
stem(mzroi_aug,points(a(n)).massSpec)

hold on 
end 
end 
linkaxes
% text(d(ind,3),d(ind,4),num2str(d(ind,5)))
%% --- Functions ---
% function clusters = referenceClustering(points,Options)
% n = length(points);
% clusters = zeros(1,n);
% refs = []; % store reference indices
% Rt_mat = [points.Rt];
% max_dist = Options.Clustering.distance_componentization;
% doPlot = Options.doPlot;
% for i = 1:n
%     assigned = false;
%     for c = 1:length(refs)
%         refIdx = refs(c);
% 
%         y = Rt_mat(:,i);%points(i).Rt;
%         x = Rt_mat(:,refIdx);%points(refIdx).Rt;
% %         dist= sqrt(sum((x - y).^2));
%         %         vecRef = [points(refIdx).elution1(:); points(refIdx).elution2(:); points(refIdx).massSpec(:)];
%         %         vecPt  = [points(i).elution1(:); points(i).elution2(:); points(i).massSpec(:)];
%         %         dist = norm(vecPt - vecRef);
%         if abs(y-x)<max_dist
%             profile_ref{1} = points(refIdx).elution1;
%             profile_cmp{1} = points(i).elution1;
%             profile_ref{2} = points(refIdx).elution2;
%             profile_cmp{2} = points(i).elution2;
% %             clf
%             [elu_profile, shifts] = align_elution_profiles(profile_ref, profile_cmp, doPlot,Options);
% 
%             %             cos1 = cosineSim(points(i).elution1, points(refIdx).elution1);
%             %             cos2 = cosineSim(points(i).elution2, points(refIdx).elution2);
%             cos1 = cosineSim(elu_profile{1},   points(refIdx).elution1);
%             cos2 = cosineSim(elu_profile{2},   points(refIdx).elution2);
%             cos3 = cosineSim(points(i).massSpec, points(refIdx).massSpec);
%             if cos1 * cos2 * cos3 > (1-Options.Clustering.cutoff)
%                 clusters(i) = c;
%                 assigned = true;
%                 break;
%             end
%         end
%     end
%     if ~assigned
%         refs(end+1) = i; % new cluster reference
%         clusters(i) = length(refs);
%     end
% end
% end
% 
% function c = cosineSim(v1,v2)
% c = dot(v1,v2)/(norm(v1)*norm(v2)+eps);
% end
% 
% function points = padMassSpec(points)
% maxLen = max(arrayfun(@(p) length(p.massSpec), points));
% for i = 1:length(points)
%     v = points(i).massSpec(:)';
%     if length(v) < maxLen
%         v = [v zeros(1,maxLen-length(v))];
%     end
%     points(i).massSpec = v;
% end
% end
