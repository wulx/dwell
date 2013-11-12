%PROFOPT motion profile optimization to obtain the minimum RMSE of dwell time
clear all; close all; clc
% add lsramp path
if isempty(strfind(path, 'lsramp'))
    folderExist = exist('lsramp', 'dir');
    if folderExist == 7 % folder
        oldpath = addpath(fullfile(pwd, 'lsramp'));
    end
end



%% optimization process

% #1 load data and pre-process -------------------------------------------%
load dwell_time.mat
% select the 3rd stroke as benchmark test
ogee = downDwellTime(:, 3);

% spline interpolation (or extrapolation) and smoothing
[xx, ys] = spis(ogee);

% low pass filtering
lowPass = ys>maxDwellTime;
if any(lowPass)
    ys(lowPass) = maxDwellTime;
end

figure, hold on;
plot(xx, ys, 'k-')

% find critical points
opts = {'MINPEAKHEIGHT', -1, 'MINPEAKDISTANCE', 15, 'NPEAKS', 5};
[locs, ptype] = findpoints(ys, opts);
crtLocs = [1 locs numel(xx)]; % all critical locations
crtTypes = [0 ptype 0]; % treate boundary points as infletion points
nCrts = numel(crtLocs);

% mark all critical points at the spline
for i = 1:nCrts
    li = crtLocs(i);
    ti = crtTypes(i);
    if ti > 0 % peaks
        plot(xx(li), ys(li), 'ro')
    elseif ti < 0 % valleys
        plot(xx(li), ys(li), 'bo')
    else % inflection points
        plot(xx(li), ys(li), 'ko')
    end
end

crtDwells = ys(crtLocs);
crtAngs = asin(crtDwells / maxDwellTime);
crtDegs = rad2deg(crtAngs);

stepAngleDeg = 1.8 / 8;
crtSteps = round(crtDegs / stepAngleDeg);

% steps difference should not be less than 4
narrows = abs(crtSteps(2:end) - crtSteps(1:end-1)) < 4;
broadenMask = [narrows(1) narrows(1:end)] | [narrows(1:end) narrows(end)];

% broaden steps difference
amp1 = 3 * crtTypes; % @param --------------------------------------------@
amp2 = 3 * broadenMask .* crtTypes; % @param -----------------------------@

crtSteps = crtSteps + amp1 + amp2;


% #2 protocol of dwell time algorithm ------------------------------------%

% set params
stroke = 460;
% pre-allocation for steps and time sequencies
S = cell(1, nCrts-1);
T = cell(1, nCrts-1);

for j = 2:numel(crtDegs)
    % total steps and total time
    sn_tot = (crtSteps(j) - crtSteps(j-1));
    
    % go backward or forward
    goBackward = false;
    if sn_tot < 0
        sn_tot = -sn_tot;
        goBackward = true;
    end
    
    s_u = 1;
    sn_tot = sn_tot/s_u;

    s_tot = xx(crtLocs(j)) - xx(crtLocs(j-1));
    t_tot = strokeTime * s_tot / stroke;
    
    % numbers of steps, [sn_a sn_c sn_d]
    w_a = rand; % @param -------------------------------------------------@
    sn_a = round(w_a * sn_tot);
    if sn_a < 1
        sn_a = 1;
    elseif sn_a > sn_tot-2
        sn_a = sn_tot - 2;
    end
        
    w_d = (1-w_a)*rand; % @param -----------------------------------------@
    sn_d = round(w_d * sn_tot);
    if sn_d < 1
        sn_d = 1;
    elseif sn_d > sn_tot-2
        sn_d = sn_tot - 2;
    end
    
    sn_c = sn_tot - sn_a - sn_d;
    
    % initial frequency f_i and maximum frequency f_m
    a_max = 60000; % benchtest result: 66333 pulse/s^2
    
    % 0 <= w_f <= 1
    w_f = rand; % @param -------------------------------------------------@
    
    % (0, (sn_tot/t_tot - 2)/(0.25*a_max*t_tot)]
    f_i = max(floor(sn_tot/t_tot - 0.25*a_max*t_tot*w_f), 2);
    
    % solve multiple equations of linear speed ramp
    f_m = sym('f_m', 'positive');
    sols = solve(2*sn_a/(f_i+f_m) + sn_c/f_m + 2*sn_d/(f_i+f_m) == t_tot, f_m);
    f_m = round( double(sols) );
    
    isFeasible = ~isempty(f_m) && (f_m >= sn_tot/t_tot) && ...
        (f_m^2-f_i^2 <= 2*a_max*min(sn_a, sn_d));
    
    if isFeasible
        sn = [sn_a, sn_c, sn_d];
        pf = round([f_i, f_m]);
        method = 'round';
        [f_list, dt_list] = time_per_step(sn, pf, s_u, method);
        
        % time sync
        while sum(dt_list(:)) > t_tot
            f_m = f_m + 1;
            [f_list, dt_list] = time_per_step(sn, [f_i, f_m], s_u, method);
        end
        
    end

    % time sequencies of every step
    timeSeqs = steptime(f_list, dt_list);
    
    %# ticks: stepper from 0
    steps = 0:numel(timeSeqs);
    
    if goBackward
        steps = crtSteps(j-1) - steps;
    else
        steps = crtSteps(j-1) + steps;
    end
    
    T{j-1} = timeSeqs;
    
    if j < nCrts
        S{j-1} = steps(1:end-1);
    else
        S{j-1} = steps;
    end
end

% concatenant
T = cell2mat(T);
S = cell2mat(S);

timeStep = 0.001; % 1 ms
totalTime = sum( T(:) );
nScan = ceil(totalTime / timeStep);

steps = zeros(1, nScan);
timeline = linspace(0, totalTime, nScan);

nStep = numel(T);
for i = 1:nStep
    steps(timeline>=sum(T(1:i-1)) & timeline<=sum(T(1:i))) = S(i);
end

leafWidth = 60;
initAngleDeg = 0;
projWidths = step2width(steps, stepAngleDeg, leafWidth, initAngleDeg);

scaleDivs = xx;
dwellTime = timecount(projWidths, strokeTime, scaleDivs);
%     figure, plot(scaleDivs, dwellTime);

filt = ~[1 isnan(ogee)' 1];
scaleDivs = scaleDivs(filt);
dwellTime = dwellTime(filt);

figure, hold on;
plot(scaleDivs, dwellTime, 'r-');
isNum = ~isnan(ogee);
plot(scaleDivs, ogee(isNum), 'k-')



%% outputs
% root-mean-square deviation of dwell time
% r = rmse(dwellTime, dwellTimeEstimate);