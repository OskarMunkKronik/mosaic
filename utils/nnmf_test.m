load('C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\NTS_feature_detections\Rts_2D.mat')
d = load(['C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles\SparseTensor\270724_09_dim_all_samples.mat']);
Options.coords.X= Rt1;
Options.coords.Y= Rt2;
mz_vec = [ 155.0896
  222.1489]%[111.0663,155.0896,222.1489];
% mz_vec = [426.4203
% 447.3958
% 448.4024
% 449.4056
% 464.3761];
clf
Z_tmp = [];
for n = 1:length(mz_vec)
[~,a] = min(abs(mzroi_aug-mz_vec(n)));
% clf
Z_tmp(:,:,n) = conv2(double(d.Z(:,:,a)),Options.Filter,'same');
plot3(Options.coords.X , Options.coords.Y, Z_tmp(:,:,n)'./max(Z_tmp(:,:,n),[],'all'),'Color',col(n,:))
hold on; % Keep the current plot
end
xL = [5 6.5];
yL = [8 10];
% xL = [43 46]
% yL = [24 27]
ind_1D = Options.coords.X(1,:) > xL(1) & Options.coords.X(1,:) < xL(2);
ind_2D = Options.coords.Y(:,1) > yL(1) & Options.coords.Y(:,1) < yL(2);
xlim(xL)
ylim(yL)
Z_small = Z_tmp(ind_1D,ind_2D,:);
clc
for n= 1:length(mz_vec)
[r(n),c(n)] = find( Z_small(:,:,n) == max(Z_small(:,:,n),[],'all'))
end 
% [r,c] = find( Z_small == max(Z_small(:,:,2),[],'all'))
%%
Options_tmp = Options;
% Options_tmp.curve_resolution

     [W, H, cg, exp_var, ~] = perform_nnmf_unfold(Z_small, Options);
     cg
   
%%
   % Perform nnmf on peaks where the peak apex is not identical 
        Z_tmp = get_elution_profiles_cluster(Z_feature, peak_borders, clusters, points, k);
        flag = find(~cellfun(@isempty, Z_tmp));
        for n = 1:length(flag)
        [~, H, cg, exp_var, ~] = perform_nnmf_unfold(Z_tmp{flag(n)}, Options);
        cluster_members = find(clusters{k} == flag(n));
        clusters{k}(cluster_members(cg>1) ) = max(clusters{k}) + 1;
        end 
