Z{1} = load("D:\people\mht541\rawData\nadine\SparseTensor\~20250923-004-nPFAS+OCBs_100ppb_splitless.D_Img01.cdf.mat");
Z{2} = load("D:\people\mht541\rawData\nadine\SparseTensor\~20250923-006-nPFAS+OCBs_200ppb_splitless.D_Img01.cdf.mat");

%%
mzTarget =235;
clf
[~,a] = min(abs(Z{1}.mzroi_aug-mzTarget));
col = 'rb';
for k =1:2
Zplot = double(Z{k}.Z(:,:,a))';

levels = logspace(log10(min(Zplot(Zplot > 0), [], 'all')), ...
                  log10(max(Zplot, [], 'all')), 500);

plot3(Options.coords.X, Options.coords.Y, Zplot,'Color',col(k))%, ...
        % levels, 'LineColor', col(k));hold on
        hold on
        % zlim([500 4000])

end 
% zL = zlim;
% zlim([10 zL(2)])
% view([0 90])
% axis tight