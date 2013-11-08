function [lc, cnum] = findpts(y, varargin)
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

% 13 in total, 3 params for spaps and the others are extra params for
% findpeaks.
narginchk(1, 13);

% parse inputs -----------------------------------------------------------%
p = inputParser;

ny = numel(y);
checkY = @(y) isvector(y)&(ny>=3);
addRequired(p, 'y', checkY);

defaultTol = 5.e-3;
addOptional(p, 'tol', defaultTol, @isreal); %

defaultM = 3; % quintic interpolation
checkM = @(m) isscalar(m) & ismember(m, [1 2 3]);
addOptional(p, 'm', defaultM, checkM);

% accept additional parameter value inputs
p.KeepUnmatched = true;

parse(p, y, varargin{:});

% smoothing spline -------------------------------------------------------%
tol = p.Results.tol;
m = p.Results.m;

x = 1:numel(y);
[sp, y_smooth] = spaps(x, y, tol, m);

sp1 = fnder(sp, 1);

% find critical points ---------------------------------------------------%
y1 = -abs( fnval(x, sp1) );
nargin

if nargin > 3
    varargin{3:end}
    [~, locs] = findpeaks(y1, varargin{3:end});
    % find peak locations
    [~, peakLocs] = findpeaks(y_smooth, varargin{3:end});
else
    [~, locs] = findpeaks(y1);
    [~, peakLocs] = findpeaks(y_smooth);
end

% identify the inflection points at ogee
sp2 = fnder(sp1, 1);
z = fnzeros(sp2);

inflLocs = round(z(1, :));
inflLocs = inflLocs(ismember(inflLocs, locs));

% valley locations
valleyLocs = setxor([inflLocs peakLocs], locs);

% outputs ----------------------------------------------------------------%
% cumulative location numbers of all three kinds of critical points
np = numel(peakLocs);
nv = numel(valleyLocs);
ni = numel(inflLocs);
cnum = [np, np+nv, np+nv+ni];

lc = [peakLocs, valleyLocs, inflLocs];

% xc = x(lc);
% yc = y(lc);

