clc
cd H:\PhD\Articles\8_LCxLC_Vandalf_BigDataSet\Part_2_Data\Scripts\
suspect_results = readtable("C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\SuspectScreening\Results\Reports\SuspectHits_Combined_2_mol_formula.xlsx",'ReadVariableNames',true,'Range','A1:AO466');
Options.hetero_atoms = {'O','N','Cl','Br','F','S'}

suspect_results = suspect_results(suspect_results.Num_frag >= 0,:);
suspect_results(~contains(suspect_results.molecular_formula,Options.hetero_atoms) , :) = [];
suspect_results(~contains(suspect_results.molecular_formula,{'C'}) , :) = [];
suspect_results = suspect_results(suspect_results.CosineForward > 0.7,:);
suspect_results = suspect_results(suspect_results.x_W_sec_ > 0.4,:);
suspect_results(suspect_results.Remove == 1,:) = [];

suspect_results = sortrows(suspect_results,["InChiKey","x_t__min_","x_t__sec_"]);
ind_remove = false(size(suspect_results,1),1);
for n = 1:size(suspect_results,1)
    suspect_results{n,"SuspectName"} = firstAlphaOnly(suspect_results{n,"SuspectName"});
end 
%%
[u,u_first,u_ind] = unique(suspect_results{:,"InChiKey"});
for n = 1:length(u)
    % rt_dif = abs( suspect_results{u_ind == n,["x_t__min_","x_t__sec_"]}  - suspect_results{u_first(n),["x_t__min_","x_t__sec_"]});
    % all(rt_dif < Options.Rtdev.*[0.45, 0.055],2)
    if sum(u_ind == n) >1
        vec = find(u_ind == n);
        [val,a] =   max(suspect_results{vec,"Num_frag"});
        vec(a) = [];
        ind_remove(vec) = true;
    end
end

for n = 2:size(suspect_results,1)
    rt_dif = abs( suspect_results{n,["x_t__min_","x_t__sec_"]}  - suspect_results{n-1,["x_t__min_","x_t__sec_"]}) ;
    if isequal(suspect_results{n,"InChiKey"},  suspect_results{n-1,"InChiKey"} ) && all(rt_dif < Options.Rtdev.*[0.45, 0.055])

    ind_remove(n) = true;
    end
end
suspect_results(ind_remove,:) = [];
%%
% clc
measured = suspect_results{:,"Measured_M_H__"};
all_values = cellfun(@(c) str2double(strsplit(c, ',')), ...
                     suspect_results{:, "m_zFragments"}, ...
                     'UniformOutput', false);

ppm_dev = cellfun(@(vals, m) (vals - m)./m*1e6, ...
                      all_values, num2cell(measured), ...
                      'UniformOutput', false);
frag_valid = cell2mat(cellfun(@(ppm) sum(abs(ppm) > 5)>=1 , ...
                      ppm_dev, ...
                      'UniformOutput', false));
frag_valid(:,2) = cell2mat(cellfun(@(ppm) sum(abs(ppm) > 5) >= 2, ...
                      ppm_dev, ...
                      'UniformOutput', false));
suspect_results.Frag_Valid = frag_valid;
 sum(accumarray(suspect_results.duplicateFlagInCollisionEnergyAll,frag_valid(:,1))>0)
  sum(accumarray(suspect_results.duplicateFlagInCollisionEnergyAll,frag_valid(:,1))>1)

% writetable(suspect_results, ...
%     "C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\SuspectScreening\Results\Reports\SuspectHits_Combined_2_with_frag_valid.xlsx", ...
%     'WriteVariableNames', true);

%%
% suspect_hits
rt_tmp = [Options.coords.X(1,feature_table(:,3))',Options.coords.Y(feature_table(:,4),1)];
rt_dev = [0.45*2,0.5];
Options.ppm_dev

%%
med_mz = accumarray(feature_table(:,2),feature_table(:,1),[],@median); 
feature_table(:,1) = med_mz(feature_table(:,2));
[valid_hits] = (false(size(rt_tmp,1),1));
mz_det = feature_table(:,1);
suspect_matches = cell(size(suspect_results,1),1);
validSuspectMatches = [];false(size(suspect_results,1),1);
for n = 1:size(suspect_results,1)
    % if suspect_hits
    rt_target = [suspect_results.x_t__min_(n), suspect_results.x_t__sec_(n)] ;

    valid_rt = sum([rt_tmp-rt_target]<rt_dev ,2) == 2 ;

    valid_mz = abs((mz_det - suspect_results.Measured_M_H__(n))./suspect_results.Measured_M_H__(n)*1e6) < Options.ppm_dev;

    valid_hits(valid_rt & valid_mz) = true;


    % suspect_matches{n} = unique(feature_table(valid_rt & valid_mz,2));
    % % Store the valid suspect matches for further analysis

    validSuspectMatches(valid_rt & valid_mz) = n ;
end
validSuspectMatches = validSuspectMatches(valid_hits);
% end
%%
% clf
clc
% figure(1)
n = 1;
% subplot(2,2,3)
% clf
rt_dev = [0.45*2,0.5];
col_ind = repelem([clusters_05.cluster_id],(cellfun(@numel,{clusters_05.members})));
col_ind = col_ind(order);
% final_features(order([392,395]))
scatter(coeff(:,pc_select(1)), coeff(:,pc_select(2)),'filled','MarkerFaceColor','b');
hold on
% scatter(coeff(validSuspectMatches,1), coeff(validSuspectMatches,2),'xr')
t = text(coeff(:,pc_select(1)), coeff(:,pc_select(2)), num2str(final_features'),'Color','k');
clust_vec = [6,22,23,28]
col_clust = 'rmgk'
for n = 1:length(clust_vec)
feature_list =clusters_05(clust_vec(n)).names; feature_table(valid_hits,2);[980,2836,2388,2932,4057];
% feature_list = unique(feature_list);
ind_f = []; % Initialize index array for features
matches =[];
count = 1;
col_match = lines(max(col_ind));
nn_include = false(length(feature_list),1);
for  nn = 1:length(final_features)
    % subplot(2,2,3)
    %    ;
    ind_f = final_features(nn) == feature_list; % ind_f = include_feature(valid_idx)== feature_list(nn);
    if any(ind_f)
        % t(nn)= text(coeff(nn,pc_select(1)), coeff(nn,pc_select(2)), num2str(final_features(nn)'),'Color','k')
        t(nn).Color = col_clust(n); 
        t(nn).BackgroundColor = 'W';
        col_match(col_ind(nn),:);
        scatter(coeff(nn,pc_select(1)), coeff(nn,pc_select(2)),'filled','MarkerFaceColor',col_clust(n))
        t_ind(nn) = true;
        matches(count) = str2double(t(nn).String);
        % matches_compount_type =
        % t(nn).String = unique(validSuspectMatches(ind_f));
        count = count + 1;

    else

        t(nn).String = '';
        t(nn).Color = 'none';

    end
 
end


% text(coeff(ind_f,1), coeff(ind_f,2), num2str(feature_list(nn_include)),'Color','r')
xlabel(['PC1 (', num2str(round(explained(1),1)), '%)'])
ylabel(['PC2 (', num2str(round(explained(2),1)), '%)'])
% text(0.05,0.93,[figLetter(2),') PCA Loadings'],'units','normalized','FontSize',25)
hold on
end
%%

