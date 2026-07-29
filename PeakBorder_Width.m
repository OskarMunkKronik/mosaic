% function [PeakBorder,PeakWidth] = PeakBorder_Width(timeVector,intensityVector,WidthHeight,ApexTime)
% intensityVector(isnan(intensityVector)) = 0;
% [~,ApexScan] = min(abs(timeVector - ApexTime));
% IntAtApex = intensityVector(ApexScan);
% PeakBorder(1) = timeVector(ApexScan-find(intensityVector(ApexScan:-1:1)<IntAtApex/WidthHeight,1,'first')+2);
% PeakBorder(2) = timeVector(ApexScan+find(intensityVector(ApexScan:end)<IntAtApex/WidthHeight,1,'first')-1);
% PeakWidth  = diff(PeakBorder);
% end
function [PeakBorder, PeakWidth,ApexTime] = PeakBorder_Width(timeVector, intensityVector, WidthHeight, ApexTime)
%PEAKBORDER_WIDTH Estimate chromatographic peak borders and width
%   [PeakBorder, PeakWidth] = PeakBorder_Width(timeVector, intensityVector, WidthHeight, ApexTime)
%
%   Inputs:
%     timeVector      - x-axis values (same length as intensityVector)
%     intensityVector - intensity trace
%     WidthHeight     - fraction of apex intensity used for border (e.g. 2 = half-height)
%     ApexTime        - expected apex position (same scale as timeVector)
%
%   Outputs:
%     PeakBorder [t_left, t_right]
%     PeakWidth  = PeakBorder(2) - PeakBorder(1)

% Replace NaNs with zero
intensityVector(isnan(intensityVector)) = 0;

% Find apex index
[~, ApexScan] = min(abs(timeVector - ApexTime));

idx = max(1,ApexScan-1):min(ApexScan+1,length(intensityVector));
add_idx = idx-ApexScan;
[~,ApexScan_mov] = max(intensityVector(idx));

% IntAtApex = intensityVector(ApexScan);
% clf 
% scatter(ApexScan,IntAtApex,'filled','red')
% hold on
ApexScan = ApexScan + add_idx(ApexScan_mov);
IntAtApex = intensityVector(ApexScan);
% plot(intensityVector)
% hold on
% scatter(ApexScan,IntAtApex,'filled','blue')
ApexTime = timeVector(ApexScan);
% Handle flat or invalid apex
if IntAtApex <= 0
    PeakBorder = [NaN, NaN];
    PeakWidth  = NaN;
    return;
end

% Threshold for border detection
thresh = IntAtApex / WidthHeight;

% --- Left border ---
leftIdx = find(intensityVector(1:ApexScan) < thresh, 1, 'last');
if isempty(leftIdx)
    leftIdx(1) = 1;
end
change = find(sign(diff(intensityVector(ApexScan:-1:1))) == 1,1,'first');
if ~isempty(change)
    leftIdx(2) = ApexScan - change +1;

end
leftIdx = max(leftIdx);
if isempty(leftIdx)
    leftTime = timeVector(1);  % fallback to start
else
    leftTime = timeVector(leftIdx);
end

% --- Right border ---
rightIdx = find(intensityVector(ApexScan:end) < thresh, 1, 'first') + ApexScan -1;
%%
if isempty(rightIdx)
    rightIdx(1) = length(timeVector);
    
end
change = find(sign(diff(intensityVector(ApexScan:end))) == 1,1,'first');
if ~isempty(change)
    rightIdx(2) = ApexScan + change - 1;
end
rightIdx = min(rightIdx);
%%
if isempty(rightIdx)
    rightTime = timeVector(end);  % fallback to end
else
    rightTime = timeVector(rightIdx);
end

PeakBorder = [leftTime, rightTime];
PeakWidth  = PeakBorder(2) - PeakBorder(1) +1;
% xline([leftIdx,rightIdx])
% pause(1)
% Safety check
if PeakWidth <= 0
    PeakBorder = [NaN, NaN];
    PeakWidth  = NaN;
end

end
