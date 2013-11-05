function etchDepth = dwell2depth(dwellTime, strokeTime, etchRate)
%DWELL2DEPTH convert dwell time to etch depth

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/28

% etchTime = totalTime - dwellTime
etchDepth = etchRate * (strokeTime - dwellTime);
