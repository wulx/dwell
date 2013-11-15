function [r, f_list, dt_list] = dwellopt(w_a, w_d, w_f, a1, a2, Crt, strokeTime, stepAngleDeg, leafWidth, timeStep)
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
amp1 = 10 * a1 * crtTypes; % @param -------------------------------------------@
amp2 = 10 * a2 * broadenMask .* crtTypes; % @param ----------------------------@

crtSteps = round(crtSteps + amp1 + amp2);

% #2 protocol of dwell time algorithm ------------------------------------%

% pre-allocation for steps and time sequencies
S = cell(1, nCrts-1);
T = cell(1, nCrts-1);

% handles
tc = makeTimecalc;
sf = makeSolvefm;

t_diff = 0;
for j = 2:nCrts
    % total steps and total time
    sn_tot = crtSteps(j) - crtSteps(j-1);
    
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
    % f_i = max(floor(sn_tot/t_tot - 0.25*a_max*t_tot*w_f(j-1)), 2)
    f_i = round(w_f(j-1)*(sn_tot/t_tot - 2) + 2);

    %# TODO: maximum frequency may be no solutions.
    % this is not exact solution, but is almost equals to it and less than it.
    f_m = sf(sn_a, sn_c, sn_d, f_i, t_tot);

    sn = [sn_a, sn_c, sn_d];
    pf = round([f_i, f_m]);
    method = 'round';
    [f_list, dt_list] = time_per_step(sn, pf, s_u, method);
    
    % time difference
    t_diff = t_tot - sum(dt_list(:));
    
    % time sequencies of every step
    timeSeqs = steptime(f_list, dt_list);
    
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



% time calculations ------------------------------------------------------%
%     function tc = makeTimecalc
%         
%         a = sym('a', 'positive');
%         n = sym('n', 'positive');
%         t = sym('t', 'positive');
%         
%         function tSol = timecalc(f_i, f_m)
%             
%             % solve multiple equations of linear speed ramp
%             aSol = solve((f_m + 0.5*a*s_u/f_m)^2 - (f_i - 0.5*a*s_u/f_i)^2 == 2*a*n*s_u, a);
%             a1 = aSol(end);
%             tSol = solve(((f_m + 0.5*a1*s_u/f_m) - (f_i - 0.5*a1*s_u/f_i))/a1 == t, t);
%             
%         end
%         
%         tc = @timecalc;
%     end

% solve maximum frequency ------------------------------------------------%
    function sf = makeSolvefm
        fm = sym('fm', 'positive');
        
        function f = solvefm(sn_a, sn_c, sn_d, fi, t)
            f = round(double(solve(2*(sn_a+sn_d)/(fi+fm) + sn_c/fm == t, fm)));
        end
        
        sf = @solvefm;
    end
        
end

% solve maximum frequency ------------------------------------------------%
% function sf = makeSolvefm
% 
% fm = sym('fm', 'positive');
% a = sym('a', 'positive');
% n = sym('n', 'positive');
% 
%     function f = solvefm(sn_a, sn_c, sn_d, fi, t, s)
%         
%         if nargin<6, s = 1; end
%         
%         %fm0 = solve(2*sn_a/(fi+fm) + sn_c/fm + 2*sn_d/(fi+fm) == t, fm);
%         
%         % select the second acceleration as defaut value
%         aa = solve((fm + 0.5*a*s/fm)^2 - (fi - 0.5*a*s/fi)^2 == 2*a*n*s, a);
%         
%         % get acceleration and deceleration by symbolic substitution
%         aa2 = aa(2);
%         acc = subs(aa2, n, sn_a);
%         dec = subs(aa2, n, sn_d);
%         
%         fm1 = solve(((fm + 0.5*acc*s/fm) - (fi - 0.5*acc*s/fi))/acc + sn_c/fm + ((fm + 0.5*dec*s/fm) - (fi - 0.5*dec*s/fi))/dec == t);
%         
%         f = round(double(fm1));
% 
%     end
% 
% sf = @solvefm;
% end
