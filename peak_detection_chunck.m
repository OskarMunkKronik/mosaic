function [peak_borders ,sum_elution_profile ,Z_feature,allPeaks ] = peak_detection_chunck(data,a_vec,Options)
allPeaks     = [];
[peak_borders,sum_elution_profile,Z_feature]  = deal(cell(length(a_vec),1));
n_feature = 1;

sz = size(data); %size
if length(sz) < 3
    sz(3) =1;
end

%     data = smoothdata(data,2,'movmean',Options.num_smooth_points(2));
% for nMZ = 1:size(data,3)

% end
n_vec = find(max(data,[],1:2)>Options.int_thresh);
for n_count = 1:length(n_vec)%sz(3)
    n = n_vec(n_count);
    if Options.doSmoothing
        data(:,:,n) = conv2(data(:,:,n),Options.Filter,'same');
    end

    %detect peaks
    [pks, locs_y, locs_x] = peaks2(data(:,:,n), ...
        'MinPeakHeight', Options.int_thresh, ...
        'MinPeakDistance',Options.min_peak_distance);

    if ~isempty(pks)
        vol_tmp = zeros(length(locs_x),1);
        for n_peak = 1:length(locs_x)

            [locs_x(n_peak), locs_y(n_peak),peak_borders{n_feature},sum_elution_profile_tmp,Z_feature{n_feature}] = get_rt_square(locs_x(n_peak), locs_y(n_peak), data(:,:,n), Options);
            %% filter out bad peaks
            cosine_sim = zeros(1,2);
            if Options.doPreFilter

                for dim = 1:2

                    ref = sum(Options.Filter',dim);
                    cand = sum(Z_feature{n_feature}',dim);
                    % clf
                    % plot(ref./max(ref))
                    % hold on
                    % plot(cand./max(cand))
                    [~,~,~,cosine_sim(dim)]= cosine_similarity_1D_align(ref,cand );

                end


                if cosine_sim(1)*cosine_sim(2) > 1-Options.Clustering.cutoff && Options.doPreFilter

                    %%
                    vol_tmp(n_peak) = sum(Z_feature{n_feature},'all');
                    for dim =1:2
                        sum_elution_profile{n_feature}{dim} = zeros(sz(dim),1);
                        scan_vec = peak_borders{n_feature}{dim}(1):peak_borders{n_feature}{dim}(2);
                        sum_elution_profile{n_feature}{dim}(scan_vec) = sum_elution_profile_tmp{dim};

                    end

                    n_feature = n_feature+1;

                end
            else
                vol_tmp(n_peak) = sum(Z_feature{n_feature},'all');
                for dim =1:2
                    sum_elution_profile{n_feature}{dim} = zeros(sz(dim),1);
                    scan_vec = peak_borders{n_feature}{dim}(1):peak_borders{n_feature}{dim}(2);
                    sum_elution_profile{n_feature}{dim}(scan_vec) = sum_elution_profile_tmp{dim};

                end

                n_feature = n_feature+1;
            end

        end


        % Concatenate results into a matrix
        temp = [pks(:), locs_y(:), locs_x(:), repmat(a_vec(n), numel(pks), 1),vol_tmp];
        allPeaks = [allPeaks; temp];
    end
end
if ~isempty(allPeaks)
    allPeaks(allPeaks(:,5)==0,:) = [];
end

% end

n_feature = n_feature-1;
peak_borders = peak_borders(1:n_feature);
sum_elution_profile = sum_elution_profile(1:n_feature);
Z_feature = Z_feature(1:n_feature) ;

end

