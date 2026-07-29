load('H:\PhD\Thesis\Figures\FullScale\workflow\peak_detection.mat')
%%
clf
subplot(2,1,1)
contour(Options.coords.X,Options.coords.Y,data')
xlabel('^1t_R (min)')
ylabel('^2t_R (sec)')
hold on
scatter(Options.coords.X(1,locs_y),Options.coords.Y(locs_x,1),'red','filled')
run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')
m = [Options.coords.X(1,min(rt_window{1})),Options.coords.Y(min(rt_window{2}),1)];
width = [Options.coords.X(1,max(rt_window{1}))-m(1),Options.coords.Y(max(rt_window{2}),1)-m(2)]
b = [m,width]
rectangle('Position',b,'EdgeColor','k','LineWidth',2,'LineStyle','--');
text(0.03,1-0.1,figLetter(1),'FontSize',Font(3),'Units','normalized')
xlim([27 35])
ylim([5 10])
% 1D elution 
lText = {'^1t_R (scan)','^2t_R (scan)'}

for dim =1:2
subplot(2,2,2+dim)u

plot(rt_window{dim},sum_elution_profile{dim},'-o')
xlabel(lText{dim})
ylabel('Summed intensity (counts)')
text(0.05,1-0.07,figLetter(1+dim),'FontSize',Font(3),'Units','normalized')
    xline(peak_borders{dim},'--r','LineWidth',2)
    xlim([min(rt_window{dim})-1 , max(rt_window{dim})+1])

run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')

end 
fig = gcf;
fig.Position= [962 300 600 700]
exportgraphics(fig,'H:\PhD\Thesis\Figures\FullScale\workflow\peak_detection_window.tif','Resolution',1000)

