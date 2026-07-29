load('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\Rts_2D.mat')



%

cd('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection')
S_Options_peak_detection_09
ProjectPBA
data_dir = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\ROI\Samples_aug\';
load([data_dir,'mzroi_aug.mat'])
Options.coords.X = Rt1;
Options.coords.Y = Rt2;
d = full([feature_groups_all.basePeak; clusters_all_samples'; feature_groups_all.Rt;feature_groups_all.sample;feature_groups_all.volume_summed]');
%%
nIS = 1;

% clf
close all

col=lines(78);
lType={'-','--'}
for nIS = 1:length(mz_IS)
    [~,a] = min(abs([d(: ,1)- mz_IS(nIS)]))
    vec = find(d(:,1) == d(a,1));
    vec
    %     clf
    for k = 1:length(vec)
        aa = vec(k);
        %     clf
        %     if sum(d(vec,5) == d(aa,5))>1

        nn_vec = find(d(vec,5) == d(aa,5));
        for nn = 1:length(nn_vec)
            aaa = vec(nn_vec(nn));
            d(aaa,5);
            %                  subplot(1,2,1)
            %         plotFeatureGroup3D(feature_groups_all,Options,aa,'color',col( d(aa,5),:),'LineStyle',lType{nn});
            %
            if nn>1
                plotFeatureGroup3D(feature_groups_all,Options,aa,'color',col(nn,:),'LineWidth',2);

                hold on
            end
            %         subplot(1,2,2)
            %
            %         plotMassSpec(feature_groups_all,mzroi_aug,aa);
            %         title(['m/z: ',num2str(mz_IS(nIS))])
            %
            %         hold on
            % %         subplot(2,2,4)
            %         stem(mzroi_aug,-double(Z(feature_groups_all(aa).Rt(1),feature_groups_all(aa).Rt(2),:)),'Marker','none','LineWidth',1,'Color','r')

        end

        %     end
    end
    %    pause
end
run('H:\PhD\Thesis\Figures\aestethics\S_aestethics.m')
grid off
fig = gcf;
% exportgraphics(fig,'H:\PhD\Thesis\Figures\FullScale\workflow\validRT.tif','Resolution',1000)

%%
% clf
t_exp  = [44.1, 25.82;..."Tris(2-butyloxyethyl)phosphate-"
    3.15, 3;...  "Melamine- 15N3"
    5.4, 9.76;..."Acetaminophen-D3"
    11.25, 12.2;...%1H-Benzotriazole-4,5,6,7-d4
    3.6, 3.58;... %METFORMIN-D6
    40.05, 17.32;..."TERBUTRYN-D5"
    3.6, 11.28;... "Gabapentin-d4"
    37.35,23.76;...."Metolachlor-d6 (propyl-d6)"
    36.9, 20.42;.... "Terbuthylazine-d5"
    44.55, 23.18;...."6PPD-quinone-D5"
    32.85, 10.28;...."Venlafaxine-d6"
    15.3, 10.22;...    "Caffeine-D9"
    31.05, 8.48;... "Metoprolol-D7"
    31.50, 18.19 % "Diuron-d6"
    ];
for nIS = 1:length(mz_IS)
    for n = 1:length(grps{nIS})
        k_vec = d(grps{nIS}(n) == d(:,2) ,5);
        for k_count = 1:length(k_vec)
            k = k_vec(k_count);
            aa = find(k == d(:,5) & grps{nIS}(n) == d(:,2));
            for aaa = 1:length(aa)

%             if n>1
                plotFeatureGroup3D(feature_groups_all,Options,aa(aaa),'color',col(n,:),'LineWidth',2);
                hold on
%             end
            end 
        end
    end
end 