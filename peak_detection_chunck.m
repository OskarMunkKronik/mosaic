% function [peak_borders ,sum_elution_profile ,Z_feature,allPeaks ] = peak_detection_chunck(data,a_vec,Options)
% allPeaks     = [];
% [peak_borders,sum_elution_profile,Z_feature]  = deal(cell(length(a_vec),1));
% n_feature = 1;
% data_raw = data;
% sz = size(data); %size
% if length(sz) < 3
%     sz(3) =1;
% end
% 
% %     data = smoothdata(data,2,'movmean',Options.num_smooth_points(2));
% % for nMZ = 1:size(data,3)
% 
% % end
% n_vec = find(max(data,[],1:2)>Options.int_thresh);
% for n_count = 1:length(n_vec)%sz(3)
%     n = n_vec(n_count);
%     if Options.doSmoothing
% 
%         data(:,:,n) = conv2(data(:,:,n),Options.Filter,'same');
%     end
% 
%     %detect peaks
%     [pks, locs_y, locs_x] = peaks2(data(:,:,n), ...
%         'MinPeakHeight', Options.int_thresh, ...
%         'MinPeakDistance',Options.min_peak_distance);
% 
%     if ~isempty(pks)
%         n_feature_start = n_feature;
%         vol_tmp = zeros(length(locs_x),1);
%         for n_peak = 1:length(locs_x)
% 
%             [locs_x(n_peak), locs_y(n_peak),peak_borders{n_feature},sum_elution_profile_tmp,Z_feature{n_feature}] = get_rt_square(locs_x(n_peak), locs_y(n_peak), data(:,:,n), Options);
%             if Options.doMFonRaw_data
%                 Z_tmp = data_raw(peak_borders{n_feature}{1}(1):peak_borders{n_feature}{1}(2),peak_borders{n_feature}{2}(1):peak_borders{n_feature}{2}(2),n);
% 
%                 sum_elution_profile_tmp{1} = sum(Z_tmp,2);
%                 sum_elution_profile_tmp{2} = sum(Z_tmp);
%             end
%             %% filter out bad peaks
%             cosine_sim = zeros(1,2);
%             if Options.doPreFilter
% 
%                 for dim = 1:2
% 
%                     ref = sum(Options.Filter',dim);
%                     cand = sum_elution_profile_tmp{dim};
%                     [~,~,~,cosine_sim(dim)]= cosine_similarity_1D_align(ref,cand );
% 
% 
%                 end
%                 % if Options.doPlot
%                 %     clf
%                 %     for dim = 1:2
%                 %
%                 %
%                 %     end
%                 % end
% 
% 
%                 if cosine_sim(1)*cosine_sim(2) > 1-Options.Clustering.cutoff && Options.doPreFilter
% 
%                     %%
%                     vol_tmp(n_peak) = sum(Z_feature{n_feature},'all');
% 
%                     for dim =1:2
%                         sum_elution_profile{n_feature}{dim} = zeros(sz(dim),1);
%                         scan_vec = peak_borders{n_feature}{dim}(1):peak_borders{n_feature}{dim}(2);
%                         sum_elution_profile{n_feature}{dim}(scan_vec) = sum_elution_profile_tmp{dim};
% 
%                     end
% 
%                     n_feature = n_feature+1;
% 
%                 end
%             else
%                 vol_tmp(n_peak) = sum(Z_feature{n_feature},'all');
%                 for dim =1:2
%                     sum_elution_profile{n_feature}{dim} = zeros(sz(dim),1);
%                     scan_vec = peak_borders{n_feature}{dim}(1):peak_borders{n_feature}{dim}(2);
%                     sum_elution_profile{n_feature}{dim}(scan_vec) = sum_elution_profile_tmp{dim};
% 
%                 end
% 
%                 n_feature = n_feature+1;
%             end
% 
%         end
%         %%
%         borders = {peak_borders{n_feature_start:n_feature-1}};
%         overlap =eye(numel(borders));
%         % col = lines(numel(borders));
% 
% 
%         for i = 1:numel(borders)-1
%             count = 0;
% 
%             for j = i+1:numel(borders)
%                  overlap(i,j) = all([sum(sum_elution_profile{n_feature_start+i-1}{1}>0 & sum_elution_profile{n_feature_start+j-1}{1}>0)>1, sum(sum_elution_profile{n_feature_start+i-1}{2}>0 & sum_elution_profile{n_feature_start+j-1}{2}>1)]);
%             end
% 
%             if Options.doPlot
%                     clf
%                 for j = i+1:numel(borders)
% 
%                     if overlap(i,j)
%                         range1_i = borders{i}{1}(1):borders{i}{1}(2);
%                         range1_j = borders{j}{1}(1):borders{j}{1}(2);
%                         range2_i = borders{i}{2}(1):borders{i}{2}(2);
%                         range2_j = borders{j}{2}(1):borders{j}{2}(2);
%                         for dim = 1:2
%                             [X,Y]=meshgrid(range1_i,range2_i);
%                             % clf
%                             plot3(X,Y,data_raw(range1_i,range2_i,n),'r','color',col(i,:),'LineWidth',2)
%                             [X,Y]=meshgrid(range1_j,range2_j);
%                             hold on
%                             plot3(X,Y,data_raw(range1_j,range2_j,n),'--','color',col(j,:),'LineWidth',2)
%                           end
%                     end
% 
%                 end
% 
%             end
%         end
%         if any(triu(overlap,1),'all')
%             org_feature_n = length(locs_y);
%             [merged_borders_n, Z_tmp_n, peak_apex_n, ...
%                 sum_elution_profile_tmp_n, max_intensity] = ...
%                 merge_overlapping_features(overlap,borders, data_raw, n);
%             vol_tmp = cell2mat(cellfun(@(x) sum(x,'all'), Z_tmp_n, ...
%                 'UniformOutput', false));
%             pks = max_intensity;
%             locs_x = peak_apex_n(1,:);
%             locs_y = peak_apex_n(2,:);
% 
%             range_idx = n_feature_start:n_feature_start+size(overlap,1)-1;
%             c = cell(size(overlap,1),1);
%             peak_borders(range_idx) = c;
%             Z_feature(range_idx) = c;
%             sum_elution_profile(range_idx) = c;
% 
% 
%             range_idx = n_feature_start:length(locs_y)+n_feature_start-1;
%             peak_borders(range_idx) = merged_borders_n;
%             Z_feature(range_idx) = Z_tmp_n;
%             sum_elution_profile(range_idx) = sum_elution_profile_tmp_n;
% 
%             n_feature = length(locs_y)+n_feature_start;
%         % else 
%             % n_feature = n_feature +1;
%         end
%         %%
% 
% 
%         % title(max(max(sum(overlap)),max(sum(overlap,2))));
%         % Concatenate results into a matrix
%         temp = [pks(:), locs_y(:), locs_x(:), repmat(a_vec(n), numel(pks), 1),vol_tmp];
%         allPeaks = [allPeaks; temp];
%     end
% end
% 
% if ~isempty(allPeaks)
%     ind = allPeaks(:,5)==0;
%     allPeaks(ind,:) = [];
% 
% end
% 
% % end
% 
% n_feature = n_feature-1;
% peak_borders = peak_borders(1:n_feature);
% sum_elution_profile = sum_elution_profile(1:n_feature);
% Z_feature = Z_feature(1:n_feature) ;
% 
% if sum(ind)>0
% peak_borders = peak_borders(~ind);
% sum_elution_profile = sum_elution_profile(~ind);
% Z_feature = Z_feature(~ind) ;
% end 
% end
% 
function [peak_borders, sum_elution_profile, Z_feature, allPeaks] = ...
    peak_detection_chunck(data, a_vec, Options)

% -------------------------------------------------------------
% Initialization
% -------------------------------------------------------------

allPeaks = [];

[peak_borders, sum_elution_profile, Z_feature] = ...
    deal(cell(length(a_vec),1));

n_feature = 1;

data_raw = data;

sz = size(data);
if numel(sz) < 3
    sz(3) = 1;
end

% -------------------------------------------------------------
% Find slices containing signal
% -------------------------------------------------------------

n_vec = find(max(data, [], 1:2) > Options.int_thresh);

% =============================================================
% Loop over samples / slices
% =============================================================

for n_count = 1:length(n_vec)

    n = n_vec(n_count);

    % ---------------------------------------------------------
    % Smoothing
    % ---------------------------------------------------------

    if Options.doSmoothing
        data(:,:,n) = conv2(data(:,:,n), Options.Filter, 'same');
    end

    % ---------------------------------------------------------
    % Detect peaks
    % ---------------------------------------------------------

    [pks_detected, locs_y_detected, locs_x_detected] = ...
        peaks2(data(:,:,n), ...
        'MinPeakHeight', Options.int_thresh, ...
        'MinPeakDistance', Options.min_peak_distance);

    if isempty(pks_detected)
        continue
    end

    % ---------------------------------------------------------
    % Start index of features belonging to this sample
    % ---------------------------------------------------------

    n_feature_start = n_feature;

    % ---------------------------------------------------------
    % Arrays containing ONLY accepted peaks
    % ---------------------------------------------------------

    pks_valid = [];
    locs_x_valid = [];
    locs_y_valid = [];
    vol_valid = [];

    % =========================================================
    % Process each detected peak
    % =========================================================

    for n_peak = 1:length(locs_x_detected)

        % -----------------------------------------------------
        % Generate peak region
        % -----------------------------------------------------

        [loc_x, loc_y, border_tmp, ...
            sum_elution_profile_tmp, Z_tmp_feature] = ...
            get_rt_square( ...
                locs_x_detected(n_peak), ...
                locs_y_detected(n_peak), ...
                data(:,:,n), ...
                Options);

        % -----------------------------------------------------
        % Use raw data for MF if requested
        % -----------------------------------------------------

        if Options.doMFonRaw_data

            Z_raw = data_raw( ...
                border_tmp{1}(1):border_tmp{1}(2), ...
                border_tmp{2}(1):border_tmp{2}(2), ...
                n);

            sum_elution_profile_tmp{1} = sum(Z_raw, 2);
            sum_elution_profile_tmp{2} = sum(Z_raw);
        end

        % -----------------------------------------------------
        % Pre-filter
        % -----------------------------------------------------

        accept_peak = true;

        if Options.doPreFilter

            cosine_sim = zeros(1,2);

            for dim = 1:2

                ref = sum(Options.Filter', dim);
                cand = sum_elution_profile_tmp{dim};

                [~,~,~,cosine_sim(dim)] = ...
                    cosine_similarity_1D_align(ref, cand);

            end

            accept_peak = ...
                cosine_sim(1) * cosine_sim(2) > ...
                1 - Options.Clustering.cutoff;
        end

        % -----------------------------------------------------
        % Reject peak
        % -----------------------------------------------------

        if ~accept_peak
            continue
        end

        % =====================================================
        % ACCEPTED PEAK
        % =====================================================

        % Store feature information at current global index
        peak_borders{n_feature} = border_tmp;
        Z_feature{n_feature} = Z_tmp_feature;

        % -----------------------------------------------------
        % Full-size elution profiles
        % -----------------------------------------------------

        for dim = 1:2

            sum_elution_profile{n_feature}{dim} = ...
                zeros(sz(dim),1);

            scan_vec = ...
                border_tmp{dim}(1):border_tmp{dim}(2);

            sum_elution_profile{n_feature}{dim}(scan_vec) = ...
                sum_elution_profile_tmp{dim};

        end

        % -----------------------------------------------------
        % Store ONLY accepted peak information
        % -----------------------------------------------------

        pks_valid(end+1,1) = pks_detected(n_peak);
        locs_x_valid(end+1,1) = loc_x;
        locs_y_valid(end+1,1) = loc_y;

        vol_valid(end+1,1) = sum(Z_tmp_feature,'all');

        % Move global feature index forward
        n_feature = n_feature + 1;

    end

    % ---------------------------------------------------------
    % No accepted peaks in this sample
    % ---------------------------------------------------------

    if isempty(pks_valid)
        continue
    end

    % =========================================================
    % Check overlap between accepted features
    % =========================================================

    nOriginal = numel(pks_valid);

    borders = peak_borders( ...
        n_feature_start:n_feature-1);

    overlap = false(nOriginal);

    for i = 1:nOriginal-1

        for j = i+1:nOriginal

            prof_i_1 = ...
                sum_elution_profile{n_feature_start+i-1}{1};

            prof_j_1 = ...
                sum_elution_profile{n_feature_start+j-1}{1};

            prof_i_2 = ...
                sum_elution_profile{n_feature_start+i-1}{2};

            prof_j_2 = ...
                sum_elution_profile{n_feature_start+j-1}{2};

            overlap_dim1 = ...
                sum(prof_i_1 > 0 & prof_j_1 > 0) > 1;

            overlap_dim2 = ...
                sum(prof_i_2 > 0 & prof_j_2 > 0) > 1;

            overlap(i,j) = ...
                overlap_dim1 && overlap_dim2;

        end
    end

    % =========================================================
    % Merge overlapping features
    % =========================================================

    if any(overlap,'all')

        [merged_borders_n, Z_tmp_n, peak_apex_n, ...
            sum_elution_profile_tmp_n, max_intensity] = ...
            merge_overlapping_features( ...
                overlap, borders, data_raw, n);

        % -----------------------------------------------------
        % Number of final features after merging
        % -----------------------------------------------------

        nMerged = numel(merged_borders_n);

        % -----------------------------------------------------
        % Correct global range of original features
        % -----------------------------------------------------

        original_idx = ...
            n_feature_start:n_feature_start+nOriginal-1;

        % -----------------------------------------------------
        % Clear original feature entries
        % -----------------------------------------------------

        peak_borders(original_idx) = {[]};
        sum_elution_profile(original_idx) = {[]};
        Z_feature(original_idx) = {[]};

        % -----------------------------------------------------
        % New range for merged features
        % -----------------------------------------------------

        merged_idx = ...
            n_feature_start:n_feature_start+nMerged-1;

        % -----------------------------------------------------
        % Store merged features
        % -----------------------------------------------------

        peak_borders(merged_idx) = merged_borders_n;

        Z_feature(merged_idx) = Z_tmp_n;

        sum_elution_profile(merged_idx) = ...
            sum_elution_profile_tmp_n;

        % -----------------------------------------------------
        % Replace peak information with merged peak information
        %
        % peak_apex_n is N x 2:
        % column 1 = x
        % column 2 = y
        % -----------------------------------------------------

        pks_valid = max_intensity(:);

        locs_x_valid = peak_apex_n(:,1);
        locs_y_valid = peak_apex_n(:,2);

        % -----------------------------------------------------
        % Calculate volume of merged features
        % -----------------------------------------------------

        vol_valid = cellfun( ...
            @(x) sum(x,'all'), ...
            Z_tmp_n);

        vol_valid = vol_valid(:);

        % -----------------------------------------------------
        % Update global feature counter
        % -----------------------------------------------------

        n_feature = ...
            n_feature_start + nMerged;

    end

%% Plot overlapping features

if Options.doPlot

    clf
    hold on
    col = lines(nOriginal);

    for i = 1:nOriginal-1

        for j = i+1:nOriginal

            if overlap(i,j)

                range1_i = ...
                    borders{i}{1}(1):borders{i}{1}(2);

                range1_j = ...
                    borders{j}{1}(1):borders{j}{1}(2);

                range2_i = ...
                    borders{i}{2}(1):borders{i}{2}(2);

                range2_j = ...
                    borders{j}{2}(1):borders{j}{2}(2);

                % Feature i
                [X,Y] = meshgrid(range1_i, range2_i);

                plot3( ...
                    X, Y, ...
                    data_raw(range1_i,range2_i,n), ...
                    'Color', col(i,:), ...
                    'LineWidth', 2);

                % Feature j
                [X,Y] = meshgrid(range1_j, range2_j);

                plot3( ...
                    X, Y, ...
                    data_raw(range1_j,range2_j,n), ...
                    '--', ...
                    'Color', col(j,:), ...
                    'LineWidth', 2);

            end
        end
    end

    xlabel('Dimension 1')
    ylabel('Dimension 2')
    zlabel('Intensity')

    title(sprintf( ...
        'Sample %d: %d overlapping pairs', ...
        n, sum(overlap,'all')))

    view(3)
    hold off

end
    % =========================================================
    % Add final features from this sample to allPeaks
    % =========================================================

    nFinal = numel(pks_valid);

    % Safety check
    if numel(locs_x_valid) ~= nFinal || ...
       numel(locs_y_valid) ~= nFinal || ...
       numel(vol_valid) ~= nFinal

        error(['Internal indexing error for sample %d: ' ...
               'pks = %d, locs_x = %d, locs_y = %d, volume = %d'], ...
            n, ...
            numel(pks_valid), ...
            numel(locs_x_valid), ...
            numel(locs_y_valid), ...
            numel(vol_valid));
    end

    % ---------------------------------------------------------
    % allPeaks rows correspond EXACTLY to feature cells
    % ---------------------------------------------------------

    temp = [ ...
        pks_valid(:), ...
        locs_y_valid(:), ...
        locs_x_valid(:), ...
        repmat(a_vec(n), nFinal, 1), ...
        vol_valid(:)];

    allPeaks = [allPeaks; temp];

end

% =============================================================
% Final cleanup
% =============================================================

n_feature = n_feature - 1;

peak_borders = peak_borders(1:n_feature);
sum_elution_profile = sum_elution_profile(1:n_feature);
Z_feature = Z_feature(1:n_feature);

% -------------------------------------------------------------
% Final consistency check
% -------------------------------------------------------------

nFeatures = numel(peak_borders);

if nFeatures ~= size(allPeaks,1)
    error(['Final indexing mismatch: ' ...
           'peak_borders = %d, allPeaks = %d'], ...
        nFeatures, size(allPeaks,1));
end

if numel(sum_elution_profile) ~= nFeatures
    error('sum_elution_profile has incorrect number of features.');
end

if numel(Z_feature) ~= nFeatures
    error('Z_feature has incorrect number of features.');
end

end