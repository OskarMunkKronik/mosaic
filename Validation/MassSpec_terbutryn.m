nIS = 6;
% figure
clf
for kk = 1:78
k=25
clf 
d = full([feature_groups_all.basePeak; clusters_all_samples'; feature_groups_all.Rt;feature_groups_all.sample;feature_groups_all.volume_summed]');
[~,a] = min(abs(mzroi_aug-mz_nonIS(nIS)))
n = find(mzroi_aug(a) == d(:,1)& d(:,5) == k)
stem(mzroi_aug,feature_groups_all(n).massSpec,'b', ...
  'marker','none')
hold on 
[~,a] = min(abs(mzroi_aug-mz_IS(nIS)))
n = find(mzroi_aug(a) == d(:,1)& d(:,5) == k)
stem(mzroi_aug,feature_groups_all(n).massSpec,'r', ...
  'marker','none')

k = kk
d = full([feature_groups_all.basePeak; clusters_all_samples'; feature_groups_all.Rt;feature_groups_all.sample;feature_groups_all.volume_summed]');
[~,a] = min(abs(mzroi_aug-mz_nonIS(nIS)))
n = find(mzroi_aug(a) == d(:,1)& d(:,5) == k,1)
if ~isempty(n)
stem(mzroi_aug,-feature_groups_all(n).massSpec,'b', ...
  'marker','none')
end
[~,a] = min(abs(mzroi_aug-mz_IS(nIS)))
n = find(mzroi_aug(a) == d(:,1)& d(:,5) == k,1)
if ~isempty(n)
stem(mzroi_aug,-feature_groups_all(n).massSpec,'r', ...
  'marker','none')
end 
xlim([170 250])
pause
end 

%%
[~,aa] = min(abs(mzroi_aug-191.113))
[~,a] = min(abs(mzroi_aug-mz_IS(6)));
clf
% subplot(2,1,1)
for k = 1:78
    clf

    plot(-MSroi_aug{k}(aa,:)./max(MSroi_aug{k}(a,:)))
    hold on
    plot(MSroi_aug{k}(a,:)./max(MSroi_aug{k}(a,:)))
    xlim([43476    46372])
    pause

end
