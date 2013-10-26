function [nSteps, timeSeqs] = stepcount(freqStages, timeStages)
%STEPCOUNT count number of steps

nStage = numel(freqStages); % stages of pulse frequencies
% number of steps for every stage of pulse frequencies
stepNos = round( freqStages .* timeStages );
nTotal = sum( stepNos ); % total step numbers

nSteps = 1:nTotal;

timeSeqs = zeros(1, nTotal);
for i = 1:nStage
    sn_i = stepNos(i);
    sn_a = sum( stepNos(1:i-1) );
    timeSeqs(sn_a + (1:sn_i)) = 1 / freqStages(i);
end
