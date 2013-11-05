function dwellTime = depth2dwell(depth, strokeTime, etchRate)
%DEPTH2DWELL convert etch depth to dwell time

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/29

dwellTime = strokeTime - depth/etchRate;