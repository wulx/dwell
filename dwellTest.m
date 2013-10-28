%% dwell test
close all; clear all; clc

% setup
if isempty(strfind(path, 'lsramp')) 
    folderExist = exist('lsramp', 'dir');
    if folderExist == 7 % folder
        oldpath = addpath(fullfile(pwd, 'lsramp'));
    end
end

sn_a = 200;
sn_c = 60;
sn_d = 100;
pf_i = 100;
pf_m = 1000;
s_u = 2;
method = 'round';
[f_list, dt_list] = time_per_step(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method);
%[f_list, dt_list] = steps_per_time(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method);

timeSeqs = steptime(f_list, dt_list);

timeStep = 0.00001;
[steps, timeline] = timesamp(timeSeqs, timeStep);

figure, plot(timeline, steps);

stepAngleDeg = 1.8 / 8; % 8 microsteps
leafWidth = 60; % in mm
initAngleDeg = 90;
projWidths = step2width(steps, stepAngleDeg, leafWidth, initAngleDeg);

scaleDivs = 0:300;
[dwellTime, vLeaf] = dwell(projWidths, timeline, scaleDivs);

figure, plot(scaleDivs, dwellTime);
title(['translational velocity: ', num2str(vLeaf), ' ; number of steps: ', num2str(numel(timeSeqs))]);

f = 1;
etchTime = timeline(end) - timeline(1);
etchDepth = depth(dwellTime, etchTime, f);

figure, plot(scaleDivs, etchDepth);

