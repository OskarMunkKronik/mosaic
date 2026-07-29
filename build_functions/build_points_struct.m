function points = build_points_struct(allPeaks, sum_elution_profile)
% BUILD_POINTS_STRUCT - Create points struct from clustering results
%
% Inputs:
%   sample_vec          - vector of sample indices
%   allPeaks            - cell array with peak info
%   clusters            - cell array of cluster assignments
%   sum_elution_profile - cell array of elution profiles
%   mass_spectrum       - cell array of mass spectra
%   mzroi_aug           - m/z values
%
% Output:
%   points - struct array with extracted features

points = struct();
i = 1;

for n = 1:size(allPeaks,1)
    % Only keep valid features

    points(i).Rt                = allPeaks(n,[2,3])';
    points(i).elution1          = sum_elution_profile{n}{1};
    points(i).elution2          = sum_elution_profile{n}{2};
    points(i).mzroi_ind         = allPeaks(n,4);

    % volume over all peaks in the same cluster
    points(i).volume_summed     = sum(allPeaks(n,5));


    i = i + 1;

end
end
