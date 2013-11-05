%% Generate figure and remove ticklabels
close all; clear all; clc

figure, hold on;
ezplot('sin(x)', [0, pi/2])
xlabel([])

xlim([0 pi/2])
ylim([0 1])
set(gca, 'xTick', 0:pi/8:pi/2, 'yTick', 0:0.2:1)
set(gca,'xTickLabel', [], 'yTickLabel', []) %Remove tick labels

box on;
title('$\sin\!\left(\mathrm{\alpha}\right)$', 'interpreter', 'latex', 'FontSize', 14)

%% Get tick mark positions
yTicks = get(gca, 'ytick');
xTicks = get(gca, 'xtick');

ax = axis; %Get left most x-position
HorizontalOffset = 0.05;

%% Reset the ytick labels in desired font
yTickLabels = 0:0.2:1;
for i = 1:length(yTicks)    
    
    %Create text box and set appropriate properties
    text(ax(1)-HorizontalOffset, yTicks(i), ['$' num2str(yTickLabels(i)) '$'],...
        'HorizontalAlignment','Right','interpreter', 'latex','FontSize', 12);
end

%% Reset the xtick labels in desired font
minY = min(yTicks);
verticalOffset = 0.06;

xTickLabels = {'0', '\frac{\pi}{8}', '\frac{\pi}{4}', '\frac{3\, \pi}{8}', '\frac{\pi}{2}'};

for j = 1:length(xTicks)
    %Create text box and set appropriate properties
    fontSize = 12;
    if j > 1
        fontSize = 16;
    end
    text(xTicks(j), minY-verticalOffset, ['$' xTickLabels{j} '$'],...
        'HorizontalAlignment','Right','interpreter', 'latex', 'FontSize', fontSize);
end

%% 
