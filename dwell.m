function dwellTime = dwell(projWidths, scaleDivs, elapsedTime)
%DWELL dwell time algorithm for the ion beam collimator using stepper
%motors.
%
% varargin:
%   projWidths   --  projected widths
%   scaleDivs    --  scale divisions on the grating substrate
%   elapsedTime  --  total elapsed time
%
% varargout:
%   dwellTime  --  dwell time

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/31


if any(~isfinite(projWidths) | projWidths<0) 
    error('projected width should be finite non-negtive.')
end

if any(~isfinite(scaleDivs) | scaleDivs<0)
    error('substrate scale divisions should be finite non-negtive.')
end

if ~isfinite(elapsedTime) || elapsedTime<0
    error('elapsed time should be finite non-negtive.')
end

nScan = numel(projWidths);

firstPos = scaleDivs(1);
lastPos = scaleDivs(end);
leafCenterPos = linspace(firstPos, lastPos, nScan);

counts = zeros(size(scaleDivs));
for i = 1:nScan
    pos = leafCenterPos(i);
    halfWidth = 0.5*projWidths(i);
    hits = scaleDivs>(pos-halfWidth) & scaleDivs<=(pos+halfWidth);
    counts = counts + hits;
end

timePerScan = elapsedTime / nScan;
dwellTime = timePerScan * counts;

end % funciton DWELL end -------------------------------------------------%
