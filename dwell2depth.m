function etchDepth = dwell2depth(dwellTime, etchTime, f)
%DEPTH distribution of etched depth

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/28

% TODO -- there are many factors to affect the final etch depth actually 
% and here f is taken as the linear factor representing all others.
etchDepth = f * (etchTime - dwellTime);
