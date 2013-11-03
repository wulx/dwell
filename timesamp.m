function varargout = timesamp(timeSeqs, timeStep)
%TIMESAMP time sampling

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/31

totalTime = sum( timeSeqs(:) );
nScan = ceil(totalTime / timeStep);

steps = zeros(1, nScan);
timeline = linspace(0, totalTime, nScan);

nStep = numel(timeSeqs);
for i = 1:nStep
    steps = steps + (timeline >= sum(timeSeqs(1:i-1)));
end

switch nargout
    case 0
        figure, plot(timeline, steps, 'k-');
    case 1
        varargout = {steps};
    case 2
        varargout = {steps, timeline};
    otherwise
        error('too many output arguments')
end
