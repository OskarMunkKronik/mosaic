ref_spec_dir = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\SuspectScreening\ref_mass_Spectra';
cd(ref_spec_dir)


refSpecList = dir();
refSpecList = refSpecList(4:end);

mz_diff_mat = cell(numel(refSpecList),1);  % preallocate as cell
q = progressParfor(numel(refSpecList));
precursor_names =  {'precursormz','exactmass'};
fields = {'precursormz','exactmass'};
parfor n_ref = 1:numel(refSpecList)
    send(q, k);
    refSpec_tmp = dir(fullfile(refSpecList(n_ref).folder, ...
        refSpecList(n_ref).name, '*.mat'));

    tmp_all = cell(numel(refSpec_tmp),1);
    for n = 1:numel(refSpec_tmp)
        S = load(fullfile(refSpec_tmp(n).folder, refSpec_tmp(n).name));
        S.peaks(:,2) = S.peaks(:,2)./max(S.peaks(:,2))*100;
        mz = S.peaks(S.peaks(:,2)>10,1);
        if sum(abs(mz-round(mz)))>0


            f = precursor_names{(find(isfield(S.metadata, precursor_names), 1))};
            if ~isempty(f)
                mz_precursor = str2double(S.metadata.(f));
                tmp_all{n} = mz_precursor - mz( mz_precursor>mz);
            else
                tmp_all{n} = NaN;
            end
        end 
        % tmp_all{n} = val - mz;

        % pairwise differences (lower triangle only)
        % diffs = tril(abs(mz - mz'), -1);
        % tmp_all{n} = unique(diffs(diffs ~= 0));   % vectorize + remove zeros
    end
    mz_diff_mat{n_ref} = vertcat(tmp_all{:});
end

% Optional: collapse everything into a single vec
mz_diff_mat = vertcat(mz_diff_mat{:});

%% Group
mz_diff_mat = sort(mz_diff_mat,'ascend');
grps = ones(numel(mz_diff_mat),1);
g_count = 1;
wU = (1+Options.ppm_dev*1e-6) * mz_diff_mat(1);
for n = 2:numel(mz_diff_mat)
    if  mz_diff_mat(n) <= wU
        grps(n) = g_count;
    else
        wU = (1+Options.ppm_dev*1e-6) * mz_diff_mat(n);
        g_count = g_count + 1;  % Increment group count
        grps(n) = g_count;
    end

end

prob = accumarray(grps,ones(numel(mz_diff_mat),1));
mz_prob.probability = prob./sum(prob);

mz_prob.mz_delta = accumarray(grps,mz_diff_mat,[],@mean);
plot(mz_diff,prob./sum(prob))

%%
% mz_diff = round(mz_diff_mat,3);
% x = unique(mz_diff);
% [gc,grps] = groupcounts(mz_diff);
% plot(grps,gc./sum(gc))
% yline(1e-4)
% 
% % h = histogram(round(mz_diff_mat,3),x);
% % % Normalize the histogram and plot it
% % h.Normalization = 'probability';
% xlabel('\Delta m/z');
% ylabel('Probability of \Deltam/z');
% title('Histogram of m/z Differences');
% %
% % [~,a] = min(abs(x-67.0593))
% % h.Values(a)
% %
% % xlim([0 200])
% mz_prob.probality = gc./sum(gc);
% mz_prob.mz_delta = grps;
save(fullfile('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\utils\NTS_featuredetection\mz_prob_filter', 'mz_diff_data.mat'), 'mz_prob','refSpecList');
% save
%%
% h = histogram(round(mz_diff_mat,3),x);
% % Normalize the histogram and plot it
% h.Normalization = 'probability';
% xlabel('Mass-to-charge ratio (m/z)');
% ylabel('Probability of \Deltam/z');
% title('Histogram of m/z Differences');
% 
% [~,a] = min(abs(x-67.0593))
% h.Values(a)
% 
% xlim(x(a)+[-+.2 0.200])
