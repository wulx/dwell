%function r = profopt(dwellTime)
%PROFOPT motion profile optimization to obtain the minimum RMSE of dwell time 

%% operation protocols of  dwell time algorithm
% #1 time per step or steps per time -------------------------------------%
sn_a = 80;
sn_c = 40;
sn_d = 60;
pf_i = 100;
pf_m = 1000;
s_u = 2;
method = 'round';
[f_list, dt_list] = time_per_step(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method);
%[f_list, dt_list] = steps_per_time(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method);

% #2 step time -----------------------------------------------------------%
timeSeqs = steptime(f_list, dt_list);

% #3 uniform sampling of elapsed time ------------------------------------%
timeStep = 0.00001;
[steps, timeline] = timesamp(timeSeqs, timeStep);

direction = 'down';
if strcmp(direction, 'down')
    steps = - steps;
end

figure, plot(timeline, steps);

% #4 step2width ----------------------------------------------------------%
stepAngleDeg = 1.8 / 8; % 8 microsteps
leafWidth = 60; % in mm
initAngleDeg = 90;
projWidths = step2width(steps, stepAngleDeg, leafWidth, initAngleDeg);

% #5 dwell time ----------------------------------------------------------%
scaleDivs = [0, (1:200)-0.5, 200];
strokeTime = timeline(end);
margins(1:2) = 0.5*leafWidth;
dwellTime = dwell(projWidths, strokeTime, scaleDivs, margins);

%% outputs
% root-mean-square deviation of dwell time
r = rmse(dwellTime, dwellTimeEstimate);

