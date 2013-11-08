%function r = profopt(dwellTime)
%PROFOPT motion profile optimization to obtain the minimum RMSE of dwell time 

%% load objective data


%% outputs
% root-mean-square deviation of dwell time
r = rmse(dwellTime, dwellTimeEstimate);

