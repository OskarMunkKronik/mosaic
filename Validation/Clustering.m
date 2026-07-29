cd C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection
% load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251009_mz_prob.mat")
load("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_251024_mz_prob.mat")
data_dir = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\ROI\Samples_aug\';
load([data_dir,'mzroi_aug.mat'])
S_Options_peak_detection_mz_prob
%%
% load('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\No_fragmentFilter_2DSmoothing_25926.mat')
[gc,grps] =deal(cell(14,1));
% nSamples = max(d(:,5)) ;
g_max = [];
det = [];
vol_mat = nan(78,14);
% d = double(full([Results.basePeak; cluster_curve_resolution; Results.Rt;Results.sample;Results.volume_summed;Results.volume_basePeak]'));
Volume = cell(length(mz_IS),1);
d = full([feature_groups_all.basePeak; clusters_all_samples; feature_groups_all.Rt;feature_groups_all.sample;feature_groups_all.volume_summed;feature_groups_all.volume_basePeak]');
% d(:,[3,4]) =
clf
for nIS = 1:length(mz_IS)
    subplot(2,7,nIS)
    if nIS == 1
    mzTarget = mz_IS(nIS)+Options.Adducts(2)-Options.mzHydrogen;%275.2350%268.191
    else    
    mzTarget = mz_IS(nIS);%275.2350%268.191
    end 
    ind = find(abs(d(:,1)-(mzTarget))./mzTarget.*1e6<Options.ppm_dev );%+Options.Adducts(2)-Options.mzHydrogen
    % d(ind,:)
    % clf
    % col = rand()
    [gc{nIS},grps{nIS}]=groupcounts(d(ind,2));

    [~,a] = max(gc{nIS});
    % grps{nIS}(a) == d(:,2)
    gc{nIS}(a)

    validRf = abs(d(ind,[3,4])-round(rt_exp(nIS,:)./[0.45, 0.055])) < Options.Rtdev ;
    if sum(sum( validRf)>0)>1
        d_tmp = d(ind,:);
        n_det = length( unique(d_tmp(sum(validRf,2) == 2,5)));
        det(nIS) = min(n_det,    sum( sum(validRf,2) == 2));
        g_max(nIS) =min(n_det,length(unique(d(grps{nIS}(a) == d(:,2),5))));
        vol_tmp = accumarray(d(grps{nIS}(a) == d(:,2),5),d(grps{nIS}(a) == d(:,2),6),[78,1]);

        Volume_basePeak{nIS} = accumarray(d(grps{nIS}(a) == d(:,2),5),d(grps{nIS}(a) == d(:,2),7),[78,1]);

        Volume{nIS} = vol_tmp;% [d_tmp(sum(validRf,2) == 2,5),d_tmp(sum(validRf,2) == 2,6)];
        vol_mat(unique(d(grps{nIS}(a) == d(:,2),5)),nIS) = vol_tmp(vol_tmp>0);
       
    end

    col = rand(sum(ind),3);
    scatter(d(ind,3),d(ind,4),[],d(ind,2))
    text(d(ind,3),d(ind,4),num2str(d(ind,2)))
    title(num2str(max(g_max(nIS))))
    clc
end
% n
%%
clf
D = load('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\Rts_2D.mat')
Options.coords.X = D.Rt1;
Options.coords.Y= D.Rt2;
n_det_double_check = [];
nIS = 2
col = hsv(14);
for nIS = find(det<75)
mzTarget = mz_IS(nIS);
% ind = find(abs(d(:,3)-34)<4 &abs(d(:,4)-174)<5)
% scatter(d(ind,3),d(ind,4),[],d(ind,2))
% text(d(ind,3),d(ind,4),num2str(d(ind,2)))
count =0;

i = 1;
n_valid =[];

for n = 1:(size(d,1))
%         n = ind(n_count)
        count = count+sum(abs( mzroi_aug(feature_groups_all(n).massSpec>0)-mzTarget)<0.01);
        rt_check = sum(abs(squeeze(median(rt_IS(:,nIS,1:2),'omitnan'))' - [ Options.coords.X(feature_groups_all(n).Rt(2),feature_groups_all(n).Rt(1)),Options.coords.Y(feature_groups_all(n).Rt(2),feature_groups_all(n).Rt(1))])<Options.Rtdev  )==2  ;
        if sum(abs( mzroi_aug(feature_groups_all(n).massSpec>0)-mzTarget)<0.001)>0 && rt_check
%            subplot(2,1,1)
%             contour(feature_groups_all(n).elution1*feature_groups_all(n).elution2','LineColor',col(nIS,:))
%             hold on
%             subplot(2,1,2)
%             stem(mzroi_aug,feature_groups_all(n).massSpec)
            n_valid(i) = n;
            i = i+1;
%              hold on
           
        end
       
end
 n_det_double_check(nIS) = length(unique([feature_groups_all(n_valid).sample]));
   
end

count
%%
clf

for nIS = 1:length(mz_IS)
    subplot(3,5,nIS)
    mzTarget = mz_IS(nIS);

    ind = find(abs(d(:,1)-mzTarget)<0.01 );
    [u, ia, ai] = unique(d(ind,2));
    col = lines(max(ai));   % use lines() for better contrast

    %     hold on
    for n = 1:length(ind)
        ind_1D = find(feature_groups_all(ind(n)).elution1);
        ind_2D = find(feature_groups_all(ind(n)).elution2);

        % Extract coords
        rt1 = Options.coords.X(1, feature_groups_all(ind(n)).Rt(1));
        rt2 = Options.coords.Y(feature_groups_all(ind(n)).Rt(2), 1);


        EIC = feature_groups_all(ind(n)).feature_elu_2D{1};
        [r,c] = find(EIC == max(EIC,[],'all'));
        rt1_vec = [ -[(r-1):-1:1]*Options.modTime/60+rt1,rt1,[r+1:size(EIC,1)]*Options.modTime/60+rt1];
        rt2_vec = [ -[(c-1):-1:1]*Options.coords.Y(1,1)+rt2,rt2,[c+1:size(EIC,2)]*Options.coords.Y(1,1)+rt2];

        [X_tmp, Y_tmp] = meshgrid(rt1_vec, rt2_vec);
        % Dimension check
        %         if size(EIC,1) ~= numel(yvals) || size(EIC,2) ~= numel(xvals)
        %             EIC = EIC(1:min(end,numel(yvals)), 1:min(end,numel(xvals)));
        %         end
        %
        %         if size(X_tmp,2) > 1
        %             ind_zeros = sum(EIC(2:end-1,:),1) ~= 0;
        plot3(X_tmp, ...
            Y_tmp, ...
            EIC', ...
            'Color', col(ai(n),:))
        %         end
        %     end
        hold on
    end
    %%


    xlabel('^1t_R'), ylabel('^2t_R')
    xlim(rt_exp(nIS,1)+[-2 2])
    ylim(rt_exp(nIS,2)+[-2 5])
    run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')
    title({sprintf(IS_mz{nIS,1}{:}), ...
        sprintf('Samples (n: %.0f)', g_max(nIS))},'fontsize',10)
    set(gca,'FontSize',10)
end

fig =gcf;
fig.Position =   [643        1130        1278         946];


exportgraphics(fig,'H:\PhD\Thesis\Figures\FullScale\workflow\Clustering.tif','Resolution',1000)
%%
clf
b = bar([det;g_max]');
ylabel('Detected in x samples / No.correctly grouped components')
set(gca,'XTickLabel',IS_mz{:,1})
run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')

fig =gcf;
fig.Position = ([643        1130        1278         946]);
ylim([0 75])
yline(75*0.9,'LineWidth',2,'Color','k')
legend(b,{'Detected in x samples','No.correctly grouped components'},'Location','southeast')
% exportgraphics(fig,'H:\PhD\Thesis\Figures\FullScale\workflow\Clustering_bar.tif','Resolution',1000)
mean([det]'./75*100)
%%
IS_mz(find(det-g_max),1)
sum(g_max)/sum(det)*100
%%
clf
lStyle = {'-','--',':'}
cos = 0:0.01:1;
col = 'rbk';
p = plot(cos,[cos;cos.^(1/2);cos.^(1/3)],'LineWidth',3);
for S = 1:length(p)
%     p(S).Marker = markerList{S}
p(S).LineStyle = lStyle{S}
p(S).Color = col(S)
end 
% hold on 
% plot(cos,)
% plot(cos,cos.^(1/3))
lText = {'\theta_N_=_1^1^/^1, LC-HRMS (within samples)',...
    '\theta_N_=_2^1^/^2, LC\timesLC-HRMS and LC\timesIMS-HRMS (within samples), LC-HRMS (across samples)',...
    '\theta_N_=_3^1^/^3, LC\timesLC-HRMS and LC\timesIMS-HRMS (across samples)'}
legend(p,lText,'Location','southeast')
xlabel('\theta_t_h_r_e_s_h_o_l_d')
ylabel('Required \theta in each dimension (\theta_n_=_1 = ... =  \theta_n_=_N)')
run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')
fig =gcf;
fig.Position =   [643        1130        1278         946];
exportgraphics(fig,'H:\PhD\Thesis\Figures\FullScale\workflow\cosine_threshold.tif','Resolution',1000)