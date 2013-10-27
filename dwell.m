function [dwellTime, vLeaf] = dwell(projWidths, timeline, scaleDivs)
%DWELL dwell time algorithm for the ion beam collimator using stepper
%motors.
%
% varargin:
%   projWidths  --  projected widths
%   timeline    --  timeline for projected widths
%   scaleDivs   --  scale divisions on the grating substrate
%
% varargout:
%   dwellTime  --  dwell time
%   vLeaf       --  translational velocity of scanning leaf

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/24


if any(~isfinite(projWidths) | projWidths<0) 
    error('projected width should be finite non-negtive.')
end

if any(~isfinite(timeline) | timeline<0) 
    error('timeline should be finite non-negtive.')
end

if any(~isfinite(scaleDivs) | scaleDivs<0)
    error('substrate scale divisions should be finite non-negtive.')
end

nScan = numel(projWidths);

margins = 0.5 * projWidths([1 end]);

firstPos = scaleDivs(1) - margins(1);
lastPos = scaleDivs(end) + margins(2);

leafCenterPos = linspace(firstPos, lastPos, nScan);

totalDist = lastPos - firstPos;
totalTime = timeline(end) - timeline(1);
vLeaf = totalDist / totalTime;

counts = zeros(size(scaleDivs));
for i = 1:nScan
    pos = leafCenterPos(i);
    halfWidth = 0.5*projWidths(i);
    hits = scaleDivs>(pos-halfWidth) & scaleDivs<=(pos+halfWidth);
    counts = counts + hits;
end

timePerScan = totalTime / nScan;
dwellTime = timePerScan * counts;

end % funciton DWELL end -------------------------------------------------%
