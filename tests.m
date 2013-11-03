%% Find a optimal solution for elapsed time in the process of ion etching
clear all; close all; clc

% % etch depth vs stroke ---------------------------------------------------%
% load etch_depth.mat
% 
% figure, plot(stroke, depth, 'k-');
% xlabel('stroke (nm)')
% ylabel('etch depth (nm)')
% 
% % figure, plot(stroke, depth, 'k-');
% % xlabel('stroke (nm)')
% % ylabel('etch depth (nm)')
% % axis equal
% 
% % etch time vs stroke ----------------------------------------------------%
% % assume that etch rate is 1 nm/s 
% etchRate = 1;
% etchTime = depth / etchRate;
% 
% figure, hold on;
% plot(stroke, etchTime, 'k-');
% xlabel('stroke (nm)')
% ylabel('etch time (s)')
% title(['etch rate: ' num2str(etchRate) ' nm/s'])
% 
% % find peaks, valleys and boundary points
% [pks, plocs] = findpeaks(etchTime);
% [vls, vlocs] = findpeaks(-etchTime);
% vls = -vls;
% bps = etchTime([1 end]);
% 
% plot(stroke(plocs), pks, 'r.')
% plot(stroke(vlocs), vls, 'b.')
% plot(stroke([1 end]), bps, 'g.')
% 
% % find maximal an minimal
% [maxEtchTime, maxIdx] = max(etchTime);
% [minEtchTime, minIdx] = min(etchTime);
% 
% plot(stroke(maxIdx), maxEtchTime, 'ro')
% plot(stroke(minIdx), minEtchTime, 'bo')
% 
% % averages
% meanEtchTime = mean(etchTime);
% %medianEtchTime = median(etchTime);
% %middleEtchTime = 0.5 * (maxEtchTime + minEtchTime);
% 
% plot(stroke, meanEtchTime*ones(size(stroke)), 'k-.')
% 
% % elapsed time -----------------------------------------------------------%
% leafWidth = 6*10^7; % in nm
% totalStroke = stroke(end) - stroke(1) + leafWidth; % in nm
% 
% % max tunable ratio of dwell time
% maxTunableRatio = leafWidth/totalStroke;
% 
% % equations:
% % totalTime = maxEtchTime + minDwellTime = minEtchTime + maxDwellTime;
% %
% % inequations:
% % minDwellTime > 0 ==> elapsedTime > maxEtchTime
% % maxDwellTime/elapsedTime < maxTunableRatio ==>
% %    elapsedTime - minEtchTime < maxTunableRatio*elapsedTime
% % or elapsedTime < minEtchTime / (1-maxTunableRatio)
% %
% % maxEtchTime < elapsedTime < minEtchTime/(1-maxTunableRatio)
% if maxEtchTime < minEtchTime/(1-maxTunableRatio)
%     elapsedTimeRange = [maxEtchTime, minEtchTime/(1-maxTunableRatio)];
%     elapsedTimeOpt = meanEtchTime / (1 - 0.5*maxTunableRatio);
%     if elapsedTimeOpt>elapsedTimeRange(1) && elapsedTimeOpt<elapsedTimeRange(2)
%         elapsedTime = elapsedTimeOpt;
%     end
% else
%     % over the range
% end
% 
% plot(stroke, maxEtchTime*ones(size(stroke)), 'k:');
% plot(stroke, minEtchTime/(1-maxTunableRatio)*ones(size(stroke)), 'k:');
% 
% plot(stroke, elapsedTime*ones(size(stroke)), 'r:')
% plot(stroke, ( 2*meanEtchTime-elapsedTime)*ones(size(stroke)), 'b:')
% 
% % dwell time vs stroke ---------------------------------------------------%
% % assume that total running time is 32 s
% dwellTime = depth2dwell(depth, elapsedTime, etchRate);
% figure, hold on;
% plot(stroke, dwellTime, 'k-')
% xlabel('stroke (nm)')
% ylabel('dwell time (s)')
% title(['elapsed time: ' num2str(elapsedTime) ' s'])
% 
% maxDwellTime = elapsedTime - minEtchTime;
% minDwellTime = elapsedTime - maxEtchTime;
% 
% meanDwellTime = mean(dwellTime);
% %middleDwellTime = 0.5 * (maxDwellTime + minDwellTime);
% plot(stroke, meanDwellTime*ones(size(stroke)), 'k-.')
% 
% % tunable range of dwell time
% maxTunableHeight = elapsedTime * maxTunableRatio;
% tunableRange = [0 maxTunableHeight];
% plot(stroke, maxTunableHeight*ones(size(stroke)), 'b:')
% plot(stroke, zeros(size(stroke)), 'r:')
% 
% plot(stroke(plocs), elapsedTime-pks, 'r.')
% plot(stroke(vlocs), elapsedTime-vls, 'b.')
% plot(stroke([1 end]), elapsedTime-bps, 'g.')
% 
% plot(stroke(minIdx), maxDwellTime, 'bo')
% plot(stroke(maxIdx), minDwellTime, 'ro')

%% processing

% dwell time vs stroke ---------------------------------------------------%
load dwell_time.mat

figure, hold on;
plot(stroke, dwellTime, 'k-');
xlabel('stroke (mm)')
ylabel('dwell time (s)')

% find peaks, valleys and boundary points
[pks, plocs] = findpeaks(dwellTime);
[vls, vlocs] = findpeaks(-dwellTime);
vls = -vls;
bps = dwellTime([1 end]);

plot(stroke(plocs), pks, 'r.')
plot(stroke(vlocs), vls, 'b.')
plot(stroke([1 end]), bps, 'k.')
plot([stroke(1) stroke(1)], [0 bps(1)], 'k:')
plot([stroke(end) stroke(end)], [0 bps(end)], 'k:')

leafWidth = 60;
margins(1:2) = 0.5 * leafWidth;

xhead = stroke(find(stroke<margins(1))) - margins(1);
xtail = stroke(find(stroke>stroke(end)-margins(2))) + margins(2);
x = [xhead, stroke, xtail];
y = spline(stroke, dwellTime, x);

plot(x, y, 'b:')
xlim(x([1 end]))

plot(x([1 end]), y([1 end]), 'g.')

% critical points
[sortedLocs, ~] = sort([plocs vlocs]);
pvLocs = sortedLocs + numel(xhead);
crtLocs = [1 pvLocs numel(x)];

plot(x(crtLocs), y(crtLocs), 'ro')

% select the 2nd S shaped curve
% loc2 = crtLocs(2):crtLocs(3);
% plot(x(loc2), y(loc2), 'r-', 'LineWidth', 3)

elapsedTime = 32.167508015297585; % sec.

totalStroke = stroke(end) - stroke(1) + leafWidth; % in mm
% max tunable ratio of dwell time
maxTunableRatio = leafWidth/totalStroke;

vLeaf = totalStroke / elapsedTime;

% maxDwellTime = leafWidth / vLeaf;
maxDwellTime = maxTunableRatio * elapsedTime;

angles = rad2deg(asin(y / maxDwellTime));
figure, hold on;
plot(x, angles)
ylim([0 90])
xlim(x([1 end]))

format_ticks(gca, '^{}', '^{\circ}')

% microStepAngle = 0.225;
% % nMicroSteps = floor((max(angles) - min(angles)) / microStepAngle);
% minAng = ceil(min(angles)/microStepAngle) * microStepAngle;
% maxAng = floor(max(angles)/microStepAngle) * microStepAngle;
% microAngs = (minAng:microStepAngle:maxAng);


figure, hold on;
plot(stroke, dwellTime, 'k-');
xlabel('stroke (mm)')
ylabel('dwell time (s)')

% peak neighbors
ploc1 = plocs(1) + (-50:50);
plot(stroke(ploc1), dwellTime(ploc1), 'b-', 'LineWidth', 2)
plot(stroke(plocs(1)), pks(1), 'r.')
plot(stroke(plocs(1)), pks(1), 'ro')

% valley neighbors
vloc1 = vlocs(1) + (-30:30);
plot(stroke(vloc1), dwellTime(vloc1), 'b-', 'LineWidth', 2)
plot(stroke(vlocs(1)), vls(1), 'r.')
plot(stroke(vlocs(1)), vls(1), 'ro')

figure, hold on;
ploc1 = plocs(1) + (-50:50);
plot(stroke(ploc1), dwellTime(ploc1), 'b-', 'LineWidth', 2)
plot(stroke(plocs(1)), pks(1), 'r.')
plot(stroke(plocs(1)), pks(1), 'ro')

figure, hold on;
vloc1 = vlocs(1) + (-30:30);
plot(stroke(vloc1), dwellTime(vloc1), 'b-', 'LineWidth', 2)
plot(stroke(vlocs(1)), vls(1), 'r.')
plot(stroke(vlocs(1)), vls(1), 'ro')

