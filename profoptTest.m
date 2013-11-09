%PROFOPT motion profile optimization to obtain the minimum RMSE of dwell time 
clear all; close all; clc

%% load objective data
load dwell_time.mat

for j = 1:8
    
    ogee = downDwellTime(:, j);
    
    [xx, ys] = spis(ogee);
    
    figure, hold on;
    plot(xx, ys, 'k-')
    
    [locs, ptype] = findpoints(ys);
    
    for i = 1:numel(locs)
        li = locs(i);
        ti = ptype(i);
        if ti > 0 % peaks
            plot(xx(li), ys(li), 'ro')
        elseif ti < 0 % valleys
            plot(xx(li), ys(li), 'bo')
        else % inflection points
            plot(xx(li), ys(li), 'ko')
        end
    end
end
%% outputs
% root-mean-square deviation of dwell time
% r = rmse(dwellTime, dwellTimeEstimate);