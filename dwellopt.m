function [r, f_list, dt_list] = dwellopt(w_a, w_d, w_f, a1, a2, Crt, strokeTime, stepAngleDeg, leafWidth, a_max, timeStep)
%DWELLOPT motion profile optimization to obtain the minimum RMSE of dwell time
% varargin:
%   @params to be optimized
%   w_a  --  weight factors of number of acceleration steps, (0, 1)
%   w_d  --  weight factors of number of deceleration steps, (0, 1)
%   w_f  --  weight factors of initial frequency, (0, 1)
%   a1   --  amplitude of enlarging peaks or valleys
%   a2   --  amplitude of braodening narrow forks
%   @params to set
%   ogee --  initial dwell time data
%   stroke  --  distance in mm
%   stepAngleDeg  --  step angle in degree
%   a_max  --  maximum acceleration
%   timeStep  --  time step
% varargout:
%   r       --  RMSE of dwell time
%   f_list  --  frequencies
%   dt_list --  time sequencies
%


% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/11/13


% #1 pre-processing ------------------------------------------------------%

% variable shortcuts
ogee = Crt.ogee;
scaleDivs = Crt.scaleDivs;

crtLocs = Crt.locs;
crtTypes = Crt.types;
crtSteps = Crt.steps;

nCrts = Crt.num;

% steps difference should not be less than 4
narrows = abs(crtSteps(2:end) - crtSteps(1:end-1)) < 4;
broadenMask = [narrows(1) narrows(1:end)] | [narrows(1:end) narrows(end)];

% broaden steps difference
amp1 = a1 * crtTypes; % @param -------------------------------------------@
amp2 = a2 * broadenMask .* crtTypes; % @param ----------------------------@

crtSteps = crtSteps + amp1 + amp2;

% #2 protocol of dwell time algorithm ------------------------------------%

% pre-allocation for steps and time sequencies
S = cell(1, nCrts-1);
T = cell(1, nCrts-1);

t_diff = 0;
for j = 2:nCrts
    % total steps and total time
    sn_tot = (crtSteps(j) - crtSteps(j-1));
    
    % go backward or forward
    goBackward = false;
    if sn_tot < 0
        sn_tot = -sn_tot;
        goBackward = true;
    end
    
    %# TODO: s_u can be integers greater than 1
    s_u = 1;
    sn_tot = sn_tot/s_u;
    
    s_tot = scaleDivs(crtLocs(j)) - scaleDivs(crtLocs(j-1));
    stroke = scaleDivs(end) - scaleDivs(1);
    t_tot = t_diff + strokeTime * s_tot / stroke;
    
    % numbers of steps, [sn_a sn_c sn_d]
    sn_a = round(w_a(j-1) * sn_tot);
    if sn_a < 1
        sn_a = 1;
    elseif sn_a > sn_tot-2
        sn_a = sn_tot - 2;
    end
    
    % make sure not over range
    sn_d = round((1-w_a(j-1))*w_d(j-1) * sn_tot);
    if sn_d < 1
        sn_d = 1;
    elseif sn_d > sn_tot-2
        sn_d = sn_tot - 2;
    end
    
    sn_c = sn_tot - sn_a - sn_d;
    
    % initial frequency f_i and maximum frequency f_m
    % (0, (sn_tot/t_tot - 2)/(0.25*a_max*t_tot)]
    f_i = max(floor(sn_tot/t_tot - 0.25*a_max*t_tot*w_f(j-1)), 2);
    
    %# TODO: maximum frequency may be no solutions.
    f_m = solveFm();
    sn = [sn_a, sn_c, sn_d];
    pf = round([f_i, f_m]);
    method = 'round';
    [f_list, dt_list] = time_per_step(sn, pf, s_u, method);
    
    % time sync
    while sum(dt_list(:)) > t_tot
        f_m = f_m + 1;
        [f_list, dt_list] = time_per_step(sn, [f_i, f_m], s_u, method);
    end
    % time difference
    t_diff = t_tot - sum(dt_list(:));
    
    % time sequencies of every step
    %timeSeqs = steptime(f_list, dt_list);
    timeSeqs = dt_list; % when s_u == 1
    
    %# tricks: stepper from 0
    steps = 0:numel(timeSeqs);
    
    if goBackward
        steps = crtSteps(j-1) - steps;
    else
        steps = crtSteps(j-1) + steps;
    end
    
    
    if j < nCrts
        S{j-1} = steps(1:end-1);
        T{j-1} = timeSeqs;
    else
        S{j-1} = [steps, steps(end)];
        T{j-1} = [timeSeqs, t_diff]; % time compensation
    end
end

% concatenant
T = cell2mat(T);
S = cell2mat(S);

% timeStep = 0.001; % 1 ms
nScan = ceil(strokeTime / timeStep);

stepsamp = zeros(1, nScan);
timeline = linspace(0, strokeTime, nScan);

nStep = numel(T);
for i = 1:nStep
    stepsamp(timeline>=sum(T(1:i-1)) & timeline<=sum(T(1:i))) = S(i);
end

projWidths = step2width(stepsamp, stepAngleDeg, leafWidth);

dwellTime = timecount(projWidths, strokeTime, scaleDivs);

filt = ~[1 isnan(ogee)' 1];
scaleDivs = scaleDivs(filt);
dwellTime = dwellTime(filt);

isNum = ~isnan(ogee);
ogee = ogee(isNum)';

% #3 root-mean-square deviation of dwell time
r = rmse(ogee, dwellTime);

if nargout > 1

    figure, hold on;
    plot(scaleDivs, dwellTime, 'r-');
    plot(scaleDivs, ogee, 'k-')
end

% solve maximum frequency ------------------------------------------------%
    function varargout = solveFm
        % solve multiple equations of linear speed ramp
        f_m = sym('f_m', 'positive');
        sols = solve(2*sn_a/(f_i+f_m) + sn_c/f_m + 2*sn_d/(f_i+f_m) == t_tot, f_m);
        f_m = round( double(sols) );
        
        isFeasible = ~isempty(f_m) && (f_m >= sn_tot/t_tot) && ...
            (f_m^2-f_i^2 <= 2*a_max*min(sn_a, sn_d));
        
        if isFeasible
            varargout = {f_m};
        end
    end

end
