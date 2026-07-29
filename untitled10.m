
S =6
for S = 1:size(IS_mz,1)
    clf
mzTarget = mz_IS(S)
[val,fG]=min(abs(mzroi_aug(allPeaks{k}(:,4))-mzTarget));
if val/mzTarget*1e6<10 
vizualize_feature_group(Z_feature{k},peak_borders{k},mzroi_aug,allPeaks{k},clusters{k},fG)
sgtitle(IS_mz{S,1})
pause
end 

end 
