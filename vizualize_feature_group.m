function  vizualize_feature_group(Z_feature,peak_borders,mzroi_aug,allPeaks,clusters,fG)
feature_vec =     find(clusters(fG) == clusters);
col = rand(length(feature_vec),3);
mass_spectrum = zeros(length(feature_vec),1);

for nn = 1:length(feature_vec)
    subplot(ceil(length(feature_vec)/2)+1,2,nn)
    n = feature_vec(nn);
        [X_coord,Y_coord] = meshgrid(peak_borders{n}{1}(1):peak_borders{n}{1}(2), ...
            peak_borders{n}{2}(1):peak_borders{n}{2}(2));
%     Z_recon  = sum_elution_profile{n}{1}*sum_elution_profile{n}{2}';
mass_spectrum(nn) = max( Z_feature{n},[],'all');
   [~, p{nn}] = contour(X_coord, Y_coord, Z_feature{n}');
    set(p{nn}, 'LineColor', col(nn,:));
  
 hold on
title(mzroi_aug(allPeaks(feature_vec(nn),4)))
end 
subplot(ceil(length(feature_vec)/2)+1,1,ceil(length(feature_vec)/2)+1)

% subplot(1,2,2)
stem(mzroi_aug(allPeaks(feature_vec,4)),mass_spectrum)
end