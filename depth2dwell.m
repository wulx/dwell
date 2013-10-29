function dwellTime = depth2dwell(depth, etchTime, g)
%DEPTH2DWELL convert etch depth to dwell time

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/28

dwellTime = etchTime - g*depth;