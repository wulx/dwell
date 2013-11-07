function [xc, yc] = findpts(x, y, varargin)
%FINDPTS find maxima, minima and special inflection points of a curve
% p.s. to divide the spline into S shaped curves, we should find all local 
% maxima, local minima and special inflection points. The inflection points
% refer specifically to those inflection points at ogee curves.
% algorithm:
%   [pts, locs] = findpeaks( -abs(fnval(fnder(sp), x)) );
% examples:
%   findpts(x, y, 1.e-2, 3, 'MINPEAKHEIGHT', -0.001)
%

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/11/7

% 14 in total, 4 params for spaps and the others are extra params for
% findpeaks.
narginchk(2, 14);

% parse inputs -----------------------------------------------------------%
p = inputParser;

nx = numel(x);
ny = numel(y);
checkX = @(x) isvector(x)&(nx>=3);
checkY = @(y) isvector(y)&(ny==nx);
addRequired(p, 'x', checkX);
addRequired(p, 'y', checkY);

defaultTol = 5.e-3;
addOptional(p, 'tol', defaultTol, @isreal); %

defaultM = 3; % quintic interpolation
checkM = @(m) isscalar(m) & ismember(m, [1 2 3]);
addOptional(p, 'm', defaultM, checkM);

% accept additional parameter value inputs
p.KeepUnmatched = true;

parse(p, x, y, varargin{:});

% smoothing spline -------------------------------------------------------%
tol = p.Results.tol;
m = p.Results.m;

y_smooth = spaps(x, y, tol, m);

y1 = fnval(x, fnder(y_smooth, 1));

% find critical points ---------------------------------------------------%
if nargin > 4
    [~, locs] = findpeaks(-abs(y1), varargin{3:end});
else
    [~, locs] = findpeaks(-abs(y1));
end


xc = x(locs);
yc = y(locs);

