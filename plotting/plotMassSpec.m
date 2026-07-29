function h = plotMassSpec(feature_groups_all, mzroi_aug, a,sign,norm, varargin)
%PLOTMASSSPEC  Plot a mass spectrum for a feature group.
%
%   h = plotMassSpec(feature_groups_all, mzroi_aug, a)
%   h = plotMassSpec(..., 'Color', 'r', 'LineWidth', 1.5, ...)
%
%   Inputs:
%       feature_groups_all : struct array containing feature groups
%       mzroi_aug          : m/z axis (vector)
%       a                  : index of feature group
%
%   Optional Name-Value Pairs:
%       Any additional stem() properties (Color, LineWidth, etc.)
%
%   Output:
%       h : handle to stem plot
%
%   Example:
%       plotMassSpec(feature_groups_all, mzroi_aug, 5, 'Color', 'k', 'LineWidth', 2);

    % Extract spectrum
    spec = feature_groups_all(a).massSpec;
    if norm
        spec =spec ./ max(spec);
    end 

    % Plot with stem
    h = stem(mzroi_aug, sign.*spec, 'Marker', 'none', 'LineWidth', 2, varargin{:});

    % Find nonzero indices
    nnz_ind = find(spec);

    % Set x-limits around signal region
    if ~isempty(nnz_ind)
%         xlim(mzroi_aug([min(nnz_ind), max(nnz_ind)])' + [-50 50])
  
        labels = arrayfun(@(x) sprintf('%.4f', x), full(mzroi_aug(nnz_ind)), 'UniformOutput', false);
        % text(mzroi_aug(nnz_ind), sign.*spec(nnz_ind), labels, ...
        %     'VerticalAlignment', 'bottom', ...
        %     'HorizontalAlignment', 'center', ...
        %     'FontSize', 8, 'Color', 'k');
%         
    end
    % Labels
    xlabel('m/z')
    ylabel('Intensity')

end
