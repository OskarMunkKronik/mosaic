function [row, col, peakVals] = fastPeakFind2D(M)
%FASTPEAKFIND2D Find local maxima in a 2D matrix
%   [row, col, peakVals] = fastPeakFind2D(M)
%
%   M : 2D matrix
%   row, col : coordinates of local maxima
%   peakVals : values at those maxima

    % Find regional maxima (logical mask)
    mask = imregionalmax(M);

    % Extract positions
    [row, col] = find(mask);
    peakVals = M(mask);
end
