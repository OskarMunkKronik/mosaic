k_vec = 1:length(matFile_list);
V_mh = nan(length(matFile_list),length(mz_IS));
for nIS = 1:length(mz_IS)
for k_count =1:length(k_vec)

    k = k_vec(k_count);
   
    [~,a] = min(abs([d(: ,1)- mz_IS(nIS)]));
%     [~,a_mz] = min(abs(mzroi_aug- 222.1489));
    vec = find(d(:,1) == d(a,1) & d(:,5) == k);
    if ~isempty(vec)
    aa = vec(1);
%     subplot(2,3,2)
  V_mh(k,nIS) = sum(feature_groups_all(aa).feature_elu_2D{1},'all');
    end 
end 
end 
%%

mean(std(V_mh,'omitnan')./mean(V_mh,'omitnan')*100)
% std(V_mh,'omitnan')./mean(V_mh,'omitnan')

include = [1,3:4,6:14];
mean(std(V_mh(:,include),'omitnan')./mean(V_mh(:,include),'omitnan')*100)
[t,p] = ttest([std(V_mh,'omitnan')./mean(V_mh,'omitnan')]',[std(area,'omitnan')./mean(area,'omitnan')]')
