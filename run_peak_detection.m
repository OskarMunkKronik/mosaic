
function [peak_borders ,sum_elution_profile ,Z_feature,allPeaks] = run_peak_detection(Z,mz_candidate,Options)
%% Initialize
if ~isfield(Options,'chunkSize')
    % Options.chunkSize = 10; % default
    chunkSize = 100;
end
% chunkSize = Options.chunkSize;
nChunks = ceil(numel(mz_candidate) / chunkSize);

peak_borders_cell = cell(nChunks, 1);
sum_elution_profile_cell = cell(nChunks, 1);
Z_feature_cell = cell(nChunks, 1);
allPeaks_cell = cell(nChunks, 1);
Z_chunks = cell(nChunks, 1);
mz_chunks = cell(nChunks, 1);

for chunk = 1:nChunks

    startIdx = (chunk - 1) * chunkSize + 1;
    idx_end = min(startIdx + chunkSize - 1, numel(mz_candidate));

    Z_chunks{chunk} = get_mz_slices( ...
        Z, ...
        mz_candidate(startIdx:idx_end) ...
    );

    mz_chunks{chunk} = mz_candidate(startIdx:idx_end);

end

for chunk = 1:nChunks

    %% Feature detection
    % [peak_borders_tmp, ...
    %  sum_elution_profile_tmp, ...
    %  Z_feature_tmp, ...
    %  allPeaks_tmp] = peak_detection_chunck( ...
    %     data, mz_candidate, Options);
  [peak_borders_cell{chunk}, ...
     sum_elution_profile_cell{chunk}, ...
     Z_feature_cell{chunk}, ...
     allPeaks_cell{chunk}] = peak_detection_chunck( ...
        double(Z_chunks{chunk}), mz_chunks{chunk}, Options);

    % Store results for this chunk
    % peak_borders_cell{chunk} = peak_borders_tmp;
    % sum_elution_profile_cell{chunk} = sum_elution_profile_tmp;
    % Z_feature_cell{chunk} = Z_feature_tmp;

    % if ~isempty(allPeaks_tmp)
    %     allPeaks_tmp(:,4) = allPeaks_tmp(:,4) + startIdx - 1;
    % end
    % 
    % allPeaks_cell{chunk} = allPeaks_tmp;
end

%% Concatenate results after parfor
peak_borders = vertcat(peak_borders_cell{:});
sum_elution_profile = vertcat(sum_elution_profile_cell{:});
Z_feature = vertcat(Z_feature_cell{:});
allPeaks = vertcat(allPeaks_cell{:});
% 
%  [peak_borders ,sum_elution_profile,Z_feature,allPeaks] =  deal([]);
% parfor startIdx = 1:Options.chunkSize:numel(mz_candidate)
%     idx_end = startIdx+Options.chunkSize-1;
%     data = get_mz_slices(Z,mz_candidate(startIdx:min(idx_end,numel(mz_candidate)))); %tmp data
% 
%     %% Feature detection
%     [peak_borders_tmp ,sum_elution_profile_tmp ,Z_feature_tmp,allPeaks_tmp] = peak_detection_chunck(data,mz_candidate,Options);
% 
%     peak_borders = [peak_borders;peak_borders_tmp];
%     sum_elution_profile = [sum_elution_profile;sum_elution_profile_tmp];
%     Z_feature = [Z_feature;Z_feature_tmp];
%     if ~isempty( allPeaks_tmp)
%     allPeaks_tmp(:,4) = allPeaks_tmp(:,4)+startIdx-1 ;
%     end 
%     allPeaks = [allPeaks;allPeaks_tmp];
% end

end