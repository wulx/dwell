function [xx, ys] = spis(ogee, tol, m)
%SPIS SPline Interpolation (or extrapolation) and Smoothing
%

if nargin < 3, m = 3; end
if nargin < 2, tol = 1.e-2; end

x = (1:numel(ogee))-0.5;
xx = [0 x numel(ogee)];

% filter margins
notNan = ~isnan(ogee);
ogee = ogee(notNan);
x = x(notNan);

% interpolation to envelop the boundaries
yy = spline(x, ogee, xx);

% spline smoothing
[~, ys] = spaps(xx, yy, tol, m);

