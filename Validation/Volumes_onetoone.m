%%
load('H:\PhD\Thesis\Figures\FullScale\workflow\QC_IS_data.mat','W','area','int','rt_IS')

run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')
%%
clf
xL = [1.5e4 1e7];
yL = [1.5e4 1e7];
plot(xL,yL,'Color','k','LineWidth',2)
hold on
V_std_mosaic = [];
nIS_vec = [1:14]%;
[V_std_mosaic,V_mean_mosaic,V_std_mosaic_basePeak,V_mean_mosaic_basePeak] = deal(zeros([1,size(area,2)])); % Initialize V_mean_mosaic
x = find(~isnan(area(:,1)));
for nIS = nIS_vec
    

    p(nIS) =plot(area(x,nIS),Volume{nIS}(x),'Marker',markerList{nIS},'MarkerSize',10,'MarkerFaceColor','auto','LineWidth',1);
%     p(nIS) =plot(x,Volume{nIS}(x),'Marker',markerList{nIS},'MarkerSize',10,'MarkerFaceColor','auto','LineWidth',1);
    V_std_mosaic(nIS) = std(Volume{nIS}(x));
    V_mean_mosaic(nIS) = mean(Volume{nIS}(x));
    V_std_mosaic_basePeak(nIS) = std(Volume_basePeak{nIS}(x));
    V_mean_mosaic_basePeak(nIS) = mean(Volume_basePeak{nIS}(x));
    

    hold on
    xlabel('[M+H]^+')
    ylabel('Mosaic - summed volume')
end
clc
hold on

set(gca,'XScale','log')

set(gca,'YScale','log')

legend(p(nIS_vec),IS_mz{nIS_vec,1},'Location','southeast' )
run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')
%%
figure
plot(V_std_mosaic./V_mean_mosaic*100)
mean(V_std_mosaic(nIS_vec)./V_mean_mosaic(nIS_vec)*100)
mean(std(area(x,nIS_vec))./mean(area(x,nIS_vec)))*100
% axis equal
% std(area,'omitnan')./mean(area,'omitnan')
% [t,p] = ttest([V_std_mosaic./V_mean_mosaic]',[std(area,'omitnan')./mean(area,'omitnan')]')
