clear all; close all; clc

load etch_depth.mat

figure, plot(stroke, depth, 'k-');
xlabel('stroke (nm)')
ylabel('etch depth (nm)')

% figure, plot(stroke, depth, 'k-');
% xlabel('stroke (nm)')
% ylabel('etch depth (nm)')
% axis equal

dwellTime = depth2dwell(depth, 32, 1); % assume that g = 1 nm/s
figure, plot(stroke, dwellTime, 'k-')
xlabel('stroke (nm)')
ylabel('dwell time (s)')

%dwellRange = @(stroke, leafWidth) [0, leafWidth/stroke];

leafWidth = 60; % in mm
totalStroke = stroke(end) - stroke(1);

dwellRange = [0, leafWidth/totalStroke];

%etchRate

depthPV = max(depth) - min(depth);


if dwellRange