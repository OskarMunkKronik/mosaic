mzTarget  = 152.0713285
mz_IS(14);
k = 2;
Standards = readtable('H:\PhD\Experiments\LargeSampleSet_Vandalf\Standards\POS_Mode_Agilent.xlsx','ReadRowNames',true,'Range','A1:P123');
mz_vec = Standards.mz;
rt_vec = [Standards.Rt_1D_min,Standards.Rt_2D_sec];


%%
clf
Options.doPlot = true;
n_samples = max([feature_groups_all.sample]);
[vol] = zeros(length(mz_vec), n_samples);
rt = nan(length(mz_vec), n_samples,2);
% col = lines(n_samples);
EIC = cell(n_samples,1);
for n_suspect =  1:length(mz_vec)
    if Options.doPlot
    clf
     sgtitle([Standards.CompoundName{n_suspect}, ', m/z: ',num2str(mz_vec(n_suspect))])
   
    end
    Standards.CompoundName(n_suspect)
    for k = [1:n_samples]
        mzTarget = mz_vec(n_suspect);
        n = get_cluster_members(feature_groups_all, clusters_all_samples, mzTarget, k,Options);

        if ~isempty(n)
            rt_tmp = nan(length(n),2);
            [validRt,vol_tmp] = deal(zeros(length(n),1));
            for n_rt = 1:length(n)
                rt_tmp(n_rt,:) = [Options.coords.X(1,feature_groups_all(n(n_rt)).Rt(1)),Options.coords.Y(feature_groups_all(n(n_rt)).Rt(2),1)];
                validRt(n_rt) = sum(abs(rt_tmp(n_rt,:)-rt_vec(n_suspect,:))< [1.5,0.5],2) == 2;
                vol_tmp(n_rt) = feature_groups_all(n(n_rt)).volume_summed;
            end
            % feature_groups_all.volume_summed
            [~,a_vol] = max(vol_tmp(validRt == 1));
            n = n(a_vol);
        end

        if ~isempty(n) & validRt
            vol(n_suspect,k) = feature_groups_all(n).volume_basePeak;
            rt(n_suspect,k,1:2)  = [Options.coords.X(1,feature_groups_all(n).Rt(1)),Options.coords.Y(feature_groups_all(n).Rt(2),1)];

            if Options.doPlot
                subplot(2,2,1)
                % if k 
                ind_WWTP = ismember(WWTP_type, SampleList.WWTP(k));
                [h,EIC{k}] = plotFeatureGroup3D(feature_groups_all, Options, n,'color',col(ind_WWTP,:));
                hold on
                subplot(2,2,3)
                h = plotMassSpec(feature_groups_all, mzroi_aug, n,1,false, 'color',col(ind_WWTP,:));
                hold on
            end

        end
    end
    if sum(vol(n_suspect,:))>0
        if Options.doPlot
        subplot(1,2,2)

        data   = vol(n_suspect,:);
        groups = SampleList.WWTP;

        % Sort by group
        [groupNames,~,groupIdx] = unique(groups);
        [groupsSorted, sortIdx] = sort(groups);
        dataSorted = data(sortIdx);
        groupIdxSorted = groupIdx(sortIdx);

        % colors = lines(numel(groupNames));

        hold on
        b = gobjects(numel(groupNames),1); % preallocate bar handles

        % Use consistent x positions
        x = 1:numel(dataSorted);

        for g = 1:numel(groupNames)
            idx = groupIdxSorted == g;
            b(g) = bar(find(idx),dataSorted(idx), ...
                'FaceColor', col(g,:),'EdgeColor',col(g,:));
        end
        % set(gca,'YScale','log')

        hold off
        xlabel('Sample index (sorted by group)')
        ylabel('Volume')
        title('Volumes grouped and colored by WWTP')
        legend(b, groupNames, 'Location', 'eastoutside')
        grid on
        fig = gcf;
        savefig(['H:\PhD\Articles\8_LCxLC_Vandalf_BigDataSet\Part_2_Data\Figures\Compounds_confirmed_w_standard\',Standards.CompoundName{n_suspect},'.fig'])
        exportgraphics(fig,['H:\PhD\Articles\8_LCxLC_Vandalf_BigDataSet\Part_2_Data\Figures\Compounds_confirmed_w_standard\',Standards.CompoundName{n_suspect},'.jpeg'])
        
        % pause
        end 
    end
end

%%
ind = sum(vol,2) > 0;
filteredResults = vol(ind, :);

% Create valid variable names ("1", "2", ..., "78")
varNames = cellstr({matFile_list.name});

% Row names must be a cell array of character vectors
rowName = Standards.CompoundName(ind);

% Now create the table
T = array2table(filteredResults, 'RowNames', rowName, 'VariableNames', varNames);
T = T([1:3,16,18:end],:);
% Convert numeric table to cell
C = [varNames; num2cell(table2array(T))];

% Append WWTP, Week, Day as extra rows
C(end+1,:) = SampleList.WWTP(:)';
C(end+1,:) = num2cell(SampleList.Week(:)');
C(end+1,:) = num2cell(SampleList.Day(:)');

% (Optional) back to table if needed
T2 = cell2table(C(2:end,:), 'VariableNames', C(1,:));
T = ([T2(end-2:end,:);T2(1:end-3,:)]);



writetable(T,'H:\PhD\Articles\8_LCxLC_Vandalf_BigDataSet\Part_2_Data\Tables\Targets.xlsx','WriteRowNames',true)
%%
T = readtable('H:\PhD\Articles\8_LCxLC_Vandalf_BigDataSet\Part_2_Data\Tables\Targets.xlsx','ReadRowNames',true)
% Compute correlation matrix
R = corr(T{:,:}', 'Rows', 'pairwise');

% Handle NaNs from constant or missing data
R(isnan(R)) = 0;

% Compute correlation distance (1 - |r|)
D = 1 - abs(R);

% Convert square distance matrix to condensed form for linkage
D_vec = squareform(D, 'tovector');

% Perform hierarchical clustering
Z = linkage(D_vec, 'average');

% Optimize leaf order
order = optimalleaforder(Z, D_vec);

% Reorder correlation matrix and labels
R_sorted = R(order, order);
labels_sorted = T.Properties.RowNames(order);

% Plot clustered heatmap
h = heatmap(labels_sorted, labels_sorted, R_sorted);
h.Title = 'Clustered Correlation Heatmap';
% h.XLabel = 'Variables';
% h.YLabel = 'Variables';
h.ColorLimits = [-1 1];
colormap(turbo);
