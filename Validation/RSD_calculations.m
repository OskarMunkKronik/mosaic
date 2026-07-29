RSD_vol = V_std_mosaic./V_mean_mosaic;
RSD_vol_basePeak = V_std_mosaic_basePeak./V_mean_mosaic_basePeak;
scatter (rt_IS(2,:,1), rt_IS(2,:,2),100,RSD_vol,"filled");
% colorbar
RSD_vol
%%
histogram(RSD_vol*100,'BinWidth',3)

%%
clf
RSD_raw = std(area(x,:))./(mean(area(x,:)));
scatter(RSD_raw(nIS_vec),RSD_vol(nIS_vec),'filled')
ylabel('Mosaic (RSD)')
xlabel('Target (RSD)')
hold on
plot([0,max([RSD_vol,RSD_raw],[],'all')],[0,max([RSD_vol,RSD_raw],[],'all')],'k')
lm = fitlm(RSD_raw,RSD_vol);
lm.plot

%%
clf
clear p
p(1) = scatter(RSD_raw(nIS_vec).*100,RSD_vol(nIS_vec).*100,100,'filled')

text(RSD_raw(nIS_vec).*100,RSD_vol(nIS_vec).*100,num2str(nIS_vec'),'color','k','Rotation',90,'FontSize',15)

hold on
% p(2) = scatter(RSD_raw.*100,RSD_vol_basePeak.*100,100,'r','filled')
% for nIS = nIS_vec
%     x = [RSD_raw(nIS).*100,RSD_raw(nIS)*100];
%     y = [RSD_vol(nIS).*100,RSD_vol_basePeak(nIS).*100];
%     plot(x,y,'--k')
% 
% end 
% text(RSD_raw.*100,RSD_vol_basePeak.*100,IS_mz{:,1},'color','r','Rotation',0,'FontSize',15)
hold on 
plot([0,100],[0,100],'k','LineWidth',2)
% p(2) = scatter(RSD_raw(nIS_vec).*100,RSD_vol_basePeak(nIS_vec).*100,100,'r','filled')
xlabel('RSD (%) - Raw data [M+H]')
ylabel('RSD (%) - Mosaic volume')
legend(p,{'Peak volume - All ions','Peak volume - Basepeak ion'},'Location','northwest','FontSize',25)
xlim([0 70])
grid on
axis square
 run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')