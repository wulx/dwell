%% dwell time optimization demo

clear all; close all; clc
% add lsramp path
if isempty(strfind(path, 'lsramp'))
    folderExist = exist('lsramp', 'dir');
    if folderExist == 7 % folder
        oldpath = addpath(fullfile(pwd, 'lsramp'));
    end
end

% load dwell time data
Dwell = load('dwell_time.mat');
upDwellTime = Dwell.upDwellTime;
downDwellTime = Dwell.downDwellTime;
strokeTime = Dwell.strokeTime;
maxDwellTime = Dwell.maxDwellTime;

ogee = downDwellTime(:, 3);

% struct to save data of critical points
Crt = struct;
Crt.ogee = ogee;

% spline interpolation (or extrapolation) and smoothing
[xx, ys] = spis(ogee);
Crt.scaleDivs = xx;

% low pass filtering
lowPass = ys>maxDwellTime;
if any(lowPass)
    ys(lowPass) = maxDwellTime;
end

% find critical points
opts = {'MINPEAKHEIGHT', -1, 'MINPEAKDISTANCE', 15, 'NPEAKS', 5};
[locs, ptype] = findpoints(ys, opts);
Crt.locs = [1 locs numel(xx)]; % all critical locations
Crt.types = [0 ptype 0]; % treate boundary points as infletion points

Crt.num = numel(Crt.locs);

crtDwells = ys(Crt.locs);
crtAngs = asin(crtDwells / maxDwellTime);
crtDegs = rad2deg(crtAngs);

stepAngleDeg = 1.8 / 8;
Crt.steps = round(crtDegs / stepAngleDeg);

figure, hold on;
plot(xx, ys, 'k-')

%mark all critical points at the spline
for i = 1:Crt.num
    li = Crt.locs(i);
    ti = Crt.types(i);
    if ti > 0 % peaks
        plot(xx(li), ys(li), 'ro')
    elseif ti < 0 % valleys
        plot(xx(li), ys(li), 'bo')
    else % inflection points
        plot(xx(li), ys(li), 'ko')
    end
end

leafWidth = 60; % 60 mm
a_max = 60000; % benchmark test result: 66333 pulse/s^2
timeStep = 0.001; % 1ms

nw = Crt.num-1;
% w_a = rand*ones(1, nw);
% w_d = rand*ones(1, nw);
% w_f = rand*ones(1, nw);
% a1 = 3;
% a2 = 3;

% [r, f_list, dt_list] = dwellopt(w_a, w_d, w_f, a1, a2, Crt, strokeTime, stepAngleDeg, leafWidth, a_max, timeStep);

% params = [w_a, w_d, w_f, a1, a2]
dwopt = @(params) dwellopt(params(1:nw), params(nw+(1:nw)), params(2*nw+(1:nw)), params(end-1), params(end), ...
     Crt, strokeTime, stepAngleDeg, leafWidth, timeStep);


params = [1/3*ones(1,2*nw), 0.4*ones(1,nw), 3*ones(size(Crt.types)), 3];

% [a, b, c] = dwopt(params);
lb = [zeros(1,3*nw), ones(size(Crt.types)), 1];   % Lower bound
ub = [ones(1,3*nw), 9*ones(size(Crt.types)), 9];  % Upper bound

% Open pool of MATLAB sessions for parallelcomputation
isOpen = matlabpool('size') > 0;
if ~isOpen
    nProcessers = getenv('NUMBER_OF_PROCESSORS');
    matlabpool('local', 3);
end

psopts = psoptimset('Cache', 'off', 'Vectorized','off', ...
    'UseParallel', 'always', 'CompletePoll', 'on', 'TolFun', 0.05, ...
    'TolMesh', 0.05, 'TolX', 0.05, 'PollMethod', 'GPSPositiveBasisNp1', ...
    'MaxMeshSize', 0.5, 'InitialMeshSize', 0.1, 'MeshAccelerator', 'on', ...
    'PlotFcns', {@psplotbestf,@psplotmeshsize,@psplotfuncount,@psplotbestx}, ...
    'PlotInterval', 1);
[params1, rmse1] = patternsearch(dwopt, params, [], [], [], [], lb, ub, psopts);

% if isOpen
%     matlabpool close
% end

% 
% T_STEP = 0.001;
% [rmse1, dwtime, stroke] = dwell_profile(fx, r, v, T_STEP, zfs1);
% f = fx.f; x = fx.x;
% figure, hold on;
% plot(x, f, 'r-');
% plot(stroke, dwtime, 'b-');
% axis([x(1) x(end) 0 2*r/v]);
% title(sprintf('Simulated Results ( RMSE: %f )', rmse1));
% legend('Original Plot', 'Simulated Plot');
% ylabel('Dwell time');
% xlabel('Stroke');
% hold off;
