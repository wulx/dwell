%PROFOPT motion profile optimization to obtain the minimum RMSE of dwell time 
clear all; close all; clc
%% load objective data
load dwells_test.mat

for i = 1:8
    
ogee = upDwells(:, i);

% plot(ogee)

x = (1:numel(ogee))-0.5;
xx = [0 x numel(ogee)];

yy = spline(x, ogee, xx);
[lc, cnum] = findpts(yy, 2.e-2, 3, 'MINPEAKHEIGHT', -0.1, 'MINPEAKDISTANCE', 10);

figure, hold on;
plot(xx, yy, 'k-')

peakLocs = lc(1:cnum(1));
valleyLocs = lc(cnum(1)+1:cnum(2));
inflLocs = lc(cnum(2)+1:cnum(3));
plot(x(peakLocs), ogee(peakLocs), 'ro')
plot(x(valleyLocs), ogee(valleyLocs), 'bo')
plot(x(inflLocs), ogee(inflLocs), 'ko')

xlim(xx([1 end]))

end
%% outputs
% root-mean-square deviation of dwell time
% r = rmse(dwellTime, dwellTimeEstimate);