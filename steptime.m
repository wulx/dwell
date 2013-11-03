function timeSeqs = steptime(freqStages, timeStages)
%STEPTIME step time

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/31

nStage = numel(freqStages); % stages of pulse frequencies
% number of steps for every stage of pulse frequencies
stepNos = round( freqStages .* timeStages );
nTotal = sum( stepNos ); % total step numbers

timeSeqs = zeros(1, nTotal);
for i = 1:nStage
    sn_i = stepNos(i);
    sn_a = sum( stepNos(1:i-1) ); % sum(stepNos(1:0)) == 0
    timeSeqs(sn_a + (1:sn_i)) = 1 / freqStages(i);
end
