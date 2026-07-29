
function [peak_borders ,sum_elution_profile ,Z_feature,allPeaks] = run_peak_detection(Z,mz_candidate,Options)
%% Initialize
if ~isfield(Options,'chunkSize')
    Options.chunkSize = 500; % default
end

 [peak_borders ,sum_elution_profile,Z_feature,allPeaks] =  deal([]);
for startIdx = 1:Options.chunkSize:numel(mz_candidate)
    idx_end = startIdx+Options.chunkSize-1;
    data = get_mz_slices(Z,mz_candidate(startIdx:min(idx_end,numel(mz_candidate)))); %tmp data

    %% Feature detection
    [peak_borders_tmp ,sum_elution_profile_tmp ,Z_feature_tmp,allPeaks_tmp] = peak_detection_chunck(data,mz_candidate,Options);

    peak_borders = [peak_borders;peak_borders_tmp];
    sum_elution_profile = [sum_elution_profile;sum_elution_profile_tmp];
    Z_feature = [Z_feature;Z_feature_tmp];
    if ~isempty( allPeaks_tmp)
    allPeaks_tmp(:,4) = allPeaks_tmp(:,4)+startIdx-1 ;
    end 
    allPeaks = [allPeaks;allPeaks_tmp];
end

end