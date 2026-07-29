function [h,EIC] = plotFeatureGroup3D(feature_groups_all, Options, n,varargin)
%PLOTFEATUREGROUP3D  Plot 2D elution profile of a feature group in 3D.
%
%   h = plotFeatureGroup3D(feature_groups_all, Options, n, col, ai)
%
%   Inputs:
%       feature_groups_all : struct array containing feature groups
%       Options            : struct with fields 'coords' and 'modTime'
%       n                  : index of feature group to plot
%       col                : color matrix (e.g. from lines or distinguishable_colors)
%       ai                 : group assignment vector for coloring
%
%   Output:
%       h : plot handle to the created line object
%
%   Example:
%       h = plotFeatureGroup3D(feature_groups_all, Options, 5, col, ai);

    % Extract elution indices
    ind_1D = find(feature_groups_all(n).elution1);
    ind_2D = find(feature_groups_all(n).elution2);

    % Extract central retention times
    rt1 = Options.coords.X(1, feature_groups_all(n).Rt(1));
    rt2 = Options.coords.Y(feature_groups_all(n).Rt(2), 1);

    % Extract elution matrix
    EIC = feature_groups_all(n).feature_elu_2D{1};%./max(feature_groups_all(n).feature_elu_2D{1},[],'all');%*feature_groups_all(n).volume_summed;

    % Find maximum signal position
    [r, c] = find(EIC == max(EIC, [], 'all'));

    % Construct retention time vectors
%     rt1_vec = [ -(r-1:-1:1)*Options.modTime/60 + rt1, rt1, (r+1:size(EIC,1))*Options.modTime/60 + rt1 ];
    rt1_vec = [ -(r-1:-1:1) * Options.modTime/60 + rt1, ...
             rt1, ...
             (1:(size(EIC,1)-r)) * Options.modTime/60 + rt1 ];

%     rt2_vec = [ -(c-1:-1:1)*Options.coords.Y(1,1) + rt2, rt2, (c+1:size(EIC,2))*Options.coords.Y(1,1).+1 + rt2 ];%fix vi
rt2_vec = [ -(c-1:-1:1) * Options.coords.Y(1,1) + rt2, ...
             rt2, ...
             (1:(size(EIC,2)-c)) * Options.coords.Y(1,1) + rt2 ];
    % Build coordinate grid
    [X_tmp, Y_tmp] = meshgrid(rt1_vec, rt2_vec);

    % Plot
    h = plot3(X_tmp, Y_tmp, EIC',varargin{:},'LineWidth',2);
    xlabel('^1t_R (min)')
    ylabel('^2t_R (sec)')
    zlabel('Component group intensity (a.u.)')
    grid on

end
