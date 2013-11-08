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

% 10 deg. ==> 75 deg
% 5 s
% 
t_tot = 2;

angs = [10 75];
ang2step = @(a, stepAngleDeg) round(a/stepAngleDeg);
sn_tot = ang2step(angs(2)-angs(1), 1.8/8);

sn_a = round(sn_tot/4); %
% sn_d = sn_a;
% sn_c = sn_tot - sn_a - sn_d;

f_mean = sn_tot / t_tot;

f_i = round(f_mean / 2); %

f_m = round(2*f_mean - f_i);

[f_list, dt_list] = time_per_step(sn, pf, s_u, method);

% sn = [80 40 60];
% pf = [100 1000];
% s_u = 2;
% 
% timeStep = 0.00001;
% 
% strokeDir = 'down';
% 
% stepAngleDeg = 1.8 / 8; % 8 microsteps
% leafWidth = 60; % in mm
% initAngleDeg = 90; % in degree
% 
% scaleDivs = [0, (1:200)-0.5, 200];
% 
% dwell(sn, pf, initAngleDeg, scaleDivs, strokeDir, s_u, leafWidth, stepAngleDeg, timeStep)
% 
