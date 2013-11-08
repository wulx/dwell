function dwell(sn, pf, initAngleDeg, scaleDivs, strokeDir, s_u, leafWidth, stepAngleDeg, timeStep, margins)
%DWELL dwell time algorithm for the ion beam collimator using stepper
%motors.
%
% varargin:
%   sn  --  cumulative stepper numbers, including sn_a, sn_c and sn_d.
%

% #1 time per step or steps per time -------------------------------------%
method = 'round';
[f_list, dt_list] = time_per_step(sn, pf, s_u, method);
%[f_list, dt_list] = steps_per_time(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method);

% #2 step time -----------------------------------------------------------%
timeSeqs = steptime(f_list, dt_list);

% #3 uniform sampling of elapsed time ------------------------------------%
[steps, timeline] = timesamp(timeSeqs, timeStep);

if strcmp(strokeDir, 'down')
    steps = - steps;
end

figure, plot(timeline, steps);

% #4 step2width ----------------------------------------------------------%
projWidths = step2width(steps, stepAngleDeg, leafWidth, initAngleDeg);

% #5 dwell time ----------------------------------------------------------%
strokeTime = timeline(end);
if nargin < 10
    margins(1:2) = 0.5*leafWidth;
end

dwellTime = timecount(projWidths, strokeTime, scaleDivs, margins);

figure, plot(scaleDivs, dwellTime);
