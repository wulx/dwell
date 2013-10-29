clear all; close all; clc

load etch_depth.mat

figure, plot(stroke, depth, 'k-');
xlabel('stroke (nm)')
ylabel('etch depth (nm)')

figure, plot(stroke, depth, 'k-');
xlabel('stroke (nm)')
ylabel('etch depth (nm)')
axis equal