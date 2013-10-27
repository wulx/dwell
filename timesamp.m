function [steps, timeline] = timesamp(timeSeqs, timeStep)
%TIMESAMP time sampling

totalTime = sum( timeSeqs(:) );
nScan = ceil(totalTime / timeStep);

steps = zeros(1, nScan);
timeline = linspace(0, totalTime, nScan);


nStep = numel(timeSeqs);
for i = 1:nStep
    steps = steps + (timeline >= sum(timeSeqs(1:i-1)));
end
    
