function [elu_profile, rt_mov_all] = align_elution_profiles(profile_ref, profile_cmp, doPlot,Options)
%ALIGN_ELUTION_PROFILES Aligns comparison profile to match reference
%
% Input:
%   profile_ref : {1x2} cell, reference elution profile (two dimensions)
%   profile_cmp : {1x2} cell, comparison elution profile (two dimensions)
%   doPlot      : logical, true to plot alignment (default = false)
%
% Output:
%   elu_profile : {1x2} cell, aligned version of profile_cmp
%   rt_mov_all  : [1x2] vector, shift amount applied in each dimension

if nargin < 3
    doPlot = false;
end

elu_profile = cell(1,2);
rt_mov_all = zeros(1,2);

if doPlot
    clf 
end 
for dim = 1:2
    % Find maxima positions
    [~,a]  = max(profile_ref{dim});
    [~,aa] = max(profile_cmp{dim});
    rt_mov = a - aa;           % how much to move cmp to match ref
    if rt_mov <= Options.Rtdev(dim)
    rt_mov_all(dim) = rt_mov;

    % Initialize aligned signal
    n = length(profile_cmp{dim});
    elu = zeros(n,1);

    % Shift comparison
    if rt_mov == 0
        elu = profile_cmp{dim};
    elseif rt_mov < 0
        elu(1:end+rt_mov) = profile_cmp{dim}(-rt_mov+1:end);
    else
        elu(rt_mov+1:end) = profile_cmp{dim}(1:end-rt_mov);
    end

    % Save aligned comparison profile
    elu_profile{dim} = elu;

    % --- Optional plotting ---
      else 
        elu_profile{dim} = profile_cmp{dim};
        elu = profile_cmp{dim};
  
  
    end
      if doPlot
        
        subplot(2,1,dim)
        hold on
        plot(profile_ref{dim}./max(profile_ref{dim}), 'b')
        plot(profile_cmp{dim}./max(profile_cmp{dim}), 'r')
        plot(elu./max(profile_cmp{dim}), '--r')
        hold off
        title(sprintf('Dim %d alignment (shift=%d)', dim, rt_mov))
        legend('Reference','Comparison','Comparison aligned')
    end
end
end
