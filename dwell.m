function dwellTime = dwell(rotAngles, timeSeqs, subScaleDivs, timeStep, vLeaf, leafWidth)
%DWELL dwell time algorithm for the ion beam collimator using stepper
%motors.
%
% varargin:
%   rotAngles     --  rotation angles
%   timeSeqs      --  time sequences
%   subScaleDivs  --  scale divisions on the grating substrate
%   vLeaf         --  translational velocity of scanning leaf
%
% varargout:
%   dwellTime  --  dwell time

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/24

if any(~isfinite(rotAngles)) 
    error('rotation angles should be finite.')
end

if any(~isfinite(timeSeqs) | timeSeqs<0) 
    error('time sequences should be finite non-negtive.')
end

if any(~isfinite(subScaleDivs) | subScaleDivs<0)
    error('substrate scale divisions should be finite non-negtive.')
end

% TODO: translational velocity of scanning leaf can be dynamically varied,
% but right now, which is constant and well meets task requirements
if ~isfinite(vLeaf) || vLeaf<0 % vLeaf is scalar
    error('translation velocity of substrate should be finite non-negtive.')
end

if ~isfinite(leafWidth) || leafWidth<0
    error('rotation angle should be finite non-negtive.')
end

% the scanning leaf parallels to ion beam when rotation angle is 0,
% and perpendicular when the angle is  pi/2.
projWidths = cos(rotAngles) * leafWidth;

% nProjWidths = numel(projWidths);

margins = 0.5 * projWidths([1 end]);

firstPos = subScaleDivs(1) - margins(1);
lastPos = subScaleDivs(end) + margins(2);

totalDist = subScaleDivs(end) - subScaleDivs(1) + sum(margins);
totalTime = totalDist / vLeaf;

% make sure take the start and end position or time point
nScan = ceil(totalTime / timeStep);
leafCenterPos = linspace(firstPos, lastPos, nScan);
timeline = linspace(0, totalTime, nScan);

counts = zeros(size(subScaleDivs));
for i = 1:nScan
    pos = leafCenterPos(i);
    halfWidth = 0.5*projWidths(i);
    hits = subScaleDivs>(pos-halfWidth) & subScaleDivs<=(pos+halfWidth);
    counts = counts + hits;
end

dwellTime = (totalTime / nScan) * counts;

end % funciton DWELL end -------------------------------------------------%
