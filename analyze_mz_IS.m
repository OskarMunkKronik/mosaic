function [det, g_max, Volume, vol_mat, gc, grps,cluster_group] = analyze_mz_IS(d, mz_IS,rt_exp, Options)
nIS = numel(mz_IS);
gc = cell(1, nIS);
grps = cell(1, nIS);
det = zeros(1, nIS);
g_max = zeros(1, nIS);
Volume = cell(1, nIS);
cluster_group = cell(1, nIS);
[vol_mat] = deal(zeros(max(d(:,5)), nIS));

for i = 1:nIS
    subplot(2,7,i)
    mzTarget = mz_IS(i);
    ind = find(abs(d(:,1) - mzTarget)./mzTarget*1e6 < Options.ppm_dev);

    if i == 2

        validRt = abs(d(ind,[3,4])-round(rt_exp(i,:)./[0.45, 0.055])) < [Options.Rtdev(1),Options.Rtdev(2)*2] ;
    else
        validRt = abs(d(ind,[3,4])-round(rt_exp(i,:)./[0.45, 0.055])) < Options.Rtdev ;
    end
    ind = ind(sum(validRt,2)  == 2);
    [gc{i}, grps{i}] = groupcounts(d(ind,2));
    [~, a] = max(gc{i});
    if ~isempty(a)

        % if sum(sum(validRt) > 0) > 1
        d_tmp = d(ind,:);
        % n_det = numel(unique(d_tmp(sum(validRt,2) == 2,5)));
        n_det = numel(unique(d_tmp(:,5)));
        det(i) = n_det;%min(n_det, sum(sum(validRt,2) == 2));
        g_max(i) = gc{i}(a);%min(n_det, numel(unique(d(grps{i}(a) == d(:,2),5))));
        vol_tmp = accumarray(d(grps{i}(a) == d(:,2),5), d(grps{i}(a) == d(:,2),6),[size(vol_mat,1),1]);
        Volume{i} = vol_tmp;
        vol_mat(:, i) = vol_tmp;
        % ind_again =any(grps{i}' == d(:,2),2);
        % validRt_again = sum(abs(d(:,[3,4])-round(rt_exp(i,:)./[0.45, 0.055])) < Options.Rtdev ,2) == 2;
        cluster_group{i} = d_tmp;%d(ind_again & validRt_again,:);
    end

    scatter(d(ind,3), d(ind,4), [], d(ind,2))
    text(d(ind,3), d(ind,4), num2str(d(ind,2)))
    title(num2str(g_max(i)))
end
% end
end
