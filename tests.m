%% Find a optimal solution for elapsed time in the process of ion etching
clear all; close all; clc

% etch depth vs stroke ---------------------------------------------------%
load etch_depth.mat

figure, plot(stroke, depth, 'k-');
xlabel('stroke (nm)')
ylabel('etch depth (nm)')

% figure, plot(stroke, depth, 'k-');
% xlabel('stroke (nm)')
% ylabel('etch depth (nm)')
% axis equal

% etch time vs stroke ----------------------------------------------------%
% assume that etch rate is 1 nm/s 
etchRate = 1;
etchTime = depth / etchRate;

figure, hold on;
plot(stroke, etchTime, 'k-');
xlabel('stroke (nm)')
ylabel('etch time (s)')

% find peaks, valleys and boundary points
[pks, plocs] = findpeaks(etchTime);
[vls, vlocs] = findpeaks(-etchTime);
vls = -vls;
bps = etchTime([1 end]);

plot(stroke(plocs), pks, 'r.')
plot(stroke(vlocs), vls, 'b.')
plot(stroke([1 end]), bps, 'g.')

% find maximal an minimal
[maxEtchTime, maxIdx] = max(etchTime);
[minEtchTime, minIdx] = min(etchTime);

plot(stroke(maxIdx), maxEtchTime, 'ro')
plot(stroke(minIdx), minEtchTime, 'bo')

% averages
meanEtchTime = mean(etchTime);
%medianEtchTime = median(etchTime);
%middleEtchTime = 0.5 * (maxEtchTime + minEtchTime);

plot(stroke, meanEtchTime*ones(size(stroke)), 'k-.')

% elapsed time -----------------------------------------------------------%
leafWidth = 6*10^7; % in nm
totalStroke = stroke(end) - stroke(1); % in nm

% max tunable ratio of dwell time
maxTunableRatio = leafWidth/totalStroke;

% equations:
% totalTime = maxEtchTime + minDwellTime = minEtchTime + maxDwellTime;
%
% inequations:
% minDwellTime > 0 ==> elapsedTime > maxEtchTime
% maxDwellTime/elapsedTime < maxTunableRatio ==>
%    elapsedTime - minEtchTime < maxTunableRatio*elapsedTime
% or elapsedTime < minEtchTime / (1-maxTunableRatio)
%
% maxEtchTime < elapsedTime < minEtchTime/(1-maxTunableRatio)
if maxEtchTime < minEtchTime/(1-maxTunableRatio)
    elapsedTimeRange = [maxEtchTime, minEtchTime/(1-maxTunableRatio)];
    elapsedTimeOpt = meanEtchTime / (1 - 0.5*maxTunableRatio);
    if elapsedTimeOpt>elapsedTimeRange(1) && elapsedTimeOpt<elapsedTimeRange(2)
        elapsedTime = elapsedTimeOpt;
    end
else
    % over the range
end

plot(stroke, maxEtchTime*ones(size(stroke)), 'k:');
plot(stroke, minEtchTime/(1-maxTunableRatio)*ones(size(stroke)), 'k:');

plot(stroke, elapsedTime*ones(size(stroke)), 'r:')
plot(stroke, ( 2*meanEtchTime-elapsedTime)*ones(size(stroke)), 'b:')

% dwell time vs stroke ---------------------------------------------------%
% assume that total running time is 32 s
dwellTime = depth2dwell(depth, elapsedTime, etchRate);
figure, hold on;
plot(stroke, dwellTime, 'k-')
xlabel('stroke (nm)')
ylabel('dwell time (s)')

maxDwellTime = elapsedTime - minEtchTime;
minDwellTime = elapsedTime - maxEtchTime;

meanDwellTime = mean(dwellTime);
%middleDwellTime = 0.5 * (maxDwellTime + minDwellTime);
plot(stroke, meanDwellTime*ones(size(stroke)), 'k-.')

% tunable range of dwell time
maxTunableHeight = elapsedTime * maxTunableRatio;
tunableRange = [0 maxTunableHeight];
plot(stroke, maxTunableHeight*ones(size(stroke)), 'b:')
plot(stroke, zeros(size(stroke)), 'r:')

plot(stroke(plocs), elapsedTime-pks, 'r.')
plot(stroke(vlocs), elapsedTime-vls, 'b.')
plot(stroke([1 end]), elapsedTime-bps, 'g.')

plot(stroke(minIdx), maxDwellTime, 'bo')
plot(stroke(maxIdx), minDwellTime, 'ro')

