function dwellTime = dwell(projWidths, timeSeqs, subScaleDivs, vLeaf)
%DWELL dwell time algorithm for the ion beam collimator using stepper
%motors.
%
% varargin:
%   projWidths    --  (a set of) projected width
%   timeSeqs      --  time sequences
%   subScaleDivs  --  scale divisions on the grating substrate
%   vSub          --  translational velocity of scanning leaf
%
% varargout:
%   dwellTime  --  dwell time

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/24


if any(~isfinite(projWidths) | projWidths<0) 
    error('all projected widths should be finite non-negtive.')
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


end % funciton DWELL end -------------------------------------------------%

function projWidths = proj(rotAngles, leafWidth)

if any(~isfinite(rotAngles))
    error('rotation angle should be finite.')
end

if ~isfinite(leafWidth) || leafWidth<0
    error('rotation angle should be finite non-negtive.')
end

projWidths = sin(rotAngles) * leafWidth;

end % funciton PROJ end --------------------------------------------------%
