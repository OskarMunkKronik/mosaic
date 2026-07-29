function [locs_x, locs_y, peak_borders,sum_elution_profile,Z_feature] = get_rt_square(locs_x, locs_y, data, Options)
%GET_RT_SQUARE  Estimate peak borders in both RT dimensions
%
%   [peak_borders] = get_rt_square(locs_x, locs_y, data, Options)
%
%   Inputs:
%     locs_x, locs_y  - peak apex coordinates
%     data            - 2D slice (matrix)
%     Options.Rtdev   - [dx, dy] window size around the apex
%
%   Output:
%     peak_borders    - cell array { [left_x, right_x], [left_y, right_y] }

    % Define retention time position
    rt = [locs_y,locs_x];
%     Z_slice = data;

    % Build RT windows around apex in both dimensions
    rt_window = cell(1,2);
    for dim = 1:2
        low  = max(1, rt(dim) - Options.Rtdev(dim));
        high = min(size(data,dim), rt(dim) + Options.Rtdev(dim));
        rt_window{dim} = low:high;
    end
    Z_feature = data(rt_window{1}, rt_window{2});

    % Summed elution profiles across each dimension
    sum_elution_profile{1} = sum(Z_feature, 2); % along rows
    sum_elution_profile{2} = sum(Z_feature, 1); % along cols

    % Estimate peak borders along each dimension
    peak_borders = cell(1,2);
    for dim = 1:2
        ApexTime = rt(dim);
        [peak_borders{dim}, ~,rt(dim)] = PeakBorder_Width(rt_window{dim}, ...
                                                  sum_elution_profile{dim}(:), ...
                                                  100, ApexTime);
    end
    locs_x = rt(1);
    locs_y = rt(2);
    Z_feature = data(peak_borders{1}(1):peak_borders{1}(2),peak_borders{2}(1):peak_borders{2}(2));
      % Summed elution profiles across each dimension
    sum_elution_profile{1} = sum(Z_feature, 2); % along rows
    sum_elution_profile{2} = sum(Z_feature, 1); % along cols

end

