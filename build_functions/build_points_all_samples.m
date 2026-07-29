function points = build_points_all_samples(fileList,Options,mzroi_aug)
%BUILD_POINTS Construct points struct from clustering and peak detection results
%
%   points = build_points(sample_vec, clusters, allPeaks, sum_elution_profile, ...
%                         Z_feature, mass_spectrum, mzroi_aug)
%
% INPUTS:
%   sample_vec          - indices of samples to process
%   clusters            - cell array of cluster assignments
%   allPeaks            - cell array of peak tables
%   sum_elution_profile - cell array of elution profiles
%   Z_feature           - cell array of feature matrices
%   mass_spectrum       - cell array of mass spectra
%   mzroi_aug           - vector of m/z values
%
% OUTPUT:
%   points              - struct array containing feature information

i = 1;
points = struct([]);
for k = 1:length(fileList)
    %     k = sample_vec(k_count);
    %     load([fileList(k).folder,'\',fileList(k).name,'*','.mat'])
    files = dir(fullfile(Options.Paths.save2mat,'MF' , [fileList(k).name '*.mat']));
    if ~isempty(files)
        load(fullfile(files(1).folder, files(1).name));  % loads the first match
    else
        error('No matching file found for %s*.mat', fileList(k).name);
    end

    for n = 1:length(clusters)
        if allPeaks(n,6) > 0
            points(i).Rt                = allPeaks(n,[2,3])';
            points(i).elution1          = sum_elution_profile{n}{1};
            points(i).elution2          = sum_elution_profile{n}{2};
            points(i).feature_elu_2D    = Z_feature(n);
            points(i).massSpec          = mass_spectra(:,clusters(n));
            points(i).volume_summed     = sum(allPeaks(clusters(n)==clusters,5));
            [points(i).volume_basePeak,a]                       = max(mass_spectra(:,clusters(n)));
            points(i).basePeak          = mzroi_aug(a);

            points(i).sample            = k;
            points(i).sample_feature_ID = allPeaks(n,6);
            points(i).sample_cluster_ID = clusters(n);

            if size(allPeaks,2) > 6
                points(i).IS_flag       = allPeaks(n,7);
            end

            i = i + 1;
        end
    end
end
end
