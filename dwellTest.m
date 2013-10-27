%% dwell test
close all; clear all; clc

% setup
if isempty(strfind(path, 'lsramp')) 
    folderExist = exist('lsramp', 'dir');
    if folderExist == 7 % folder
        oldpath = addpath(fullfile(pwd, 'lsramp'));
    end
end

sn_a = 30;
sn_c = 10;
sn_d = 18;
pf_i = 100;
pf_m = 1000;
s_u = 2;
method = 'round';
[f_list, dt_list] = time_per_step(sn_a, sn_c, sn_d, pf_i, pf_m, s_u, method);

timeSeqs = steptime(f_list, dt_list);



% stepAngleDeg = 1.8;
% radians = step2rad(nSteps, stepAngleDeg);
% 
% 
% subScaleDivs = 0:0.1:100;
% vLeaf = 200;
% leafWidth = 60;
% dwellTime = dwell(radians, timeSeqs, subScaleDivs, vLeaf, leafWidth);
% 
% plot(subScaleDivs, dwellTime);
