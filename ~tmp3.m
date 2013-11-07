%% ogee curve segmentation
close all; clear all; clc

load dwells.mat

ogee = downDwells(:, 3);
x = 1:numel(ogee);

figure, hold on;
plot(x, ogee, 'b-')

tol = 5.e-3;
w = ones(size(x)); % weighted params
m = 3; % for a quintic smoothing spline
[ogee_smooth, y] = spaps(x, ogee, tol, w, m);

%fnplt(ogee_smooth, 'r-')
plot(x, y, 'r-')
box on

hl = legend('original curve', 'smoothed spline');
set(hl, 'Location', 'SouthWest')

% smoothed ogee curve
figure;
fnplt(ogee_smooth, 'k-')
title('smoothed ogee')

% first derivative
d1 = fnder(ogee_smooth);

figure, hold on;
fnplt(d1, 'k-')
plot(x([1 end]), [0 0], 'k:')
title('first derivative')
box on

% second derivative
d2 = fnder(d1);

figure, hold on;
fnplt(d2, 'k-')
plot(x([1 end]), [0 0], 'k:')
title('second derivative')
box on

% find inflection points at ogees
v1 = fnval(d1, x);
v3 = -abs(v1);

figure, hold on;
plot(x, v3, 'k-')

[pts, locs] = findpeaks(v3);
plot(x(locs([2 end])), pts([2 end]), 'go')
plot(x(locs(1)), pts(1), 'ro')
title('find inflection points at ogees')
box on

figure, hold on;
plot(x, ogee, 'b-');
plot(x(locs(1)), ogee(locs(1)), 'ro')

title('a inflection point at the ogee curve')
box on

% z = fnzeros(fnder(ogee_smooth,2));
% z3 = z(1, :);
% 
% inflections = [];
% for i = 1:numel(locs)
%     for j = 1:numel(z3)
%         if abs(x(locs(i)) - z3(j)) < 1
%             inflections = [inflections x(locs(i))];
%         end
%     end
% end

