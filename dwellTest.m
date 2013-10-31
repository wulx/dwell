%% dwell test
close all; clear all; clc

% #0 ---------------------------------------------------------------------%
% setup
if isempty(strfind(path, 'lsramp'))
    folderExist = exist('lsramp', 'dir');
    if folderExist == 7 % folder
        oldpath = addpath(fullfile(pwd, 'lsramp'));
    end
end

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

figure, plot(timeline, steps);

% #4 step2width ----------------------------------------------------------%
stepAngleDeg = 1.8 / 8; % 8 microsteps
leafWidth = 60; % in mm
initAngleDeg = 90;
projWidths = step2width(steps, stepAngleDeg, leafWidth, initAngleDeg);

% #4 dwell time ----------------------------------------------------------%
scaleDivs = 60:300;
elapsedTime = timeline(end);
dwellTime = dwell(projWidths, scaleDivs, elapsedTime);

figure, plot(scaleDivs, dwellTime);
title(['translational velocity: ', num2str(vLeaf), ' ; number of steps: ', num2str(numel(timeSeqs))]);

