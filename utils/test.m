points = build_points_struct(allPeaks{k}, sum_elution_profile{k});
clusters = referenceClustering_one_sample(points,Options);