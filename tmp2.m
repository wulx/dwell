%% tansform and slice dwell-time map
close all; clear all; clc

% load etch depth data
load etch_depth2.mat
z = depth; % in nm
x = 0.5:399.5; % mm, 400 points in total
y = 0.5:399.5; % mm
figure, mesh(x, y, z);
axis equal;

% shear and divide projected area
ms = 2; % microsteps per full step
ionBeamWidth = 60; % mm
leafWidth = 60; % mm
subHeight = 400; % mm
subLength = 400; % mm
vStroke = leafWidth + subHeight; % vertical stroke

shx = (ionBeamWidth/ms) / vStroke;
% shx = 0.025;

nUpStroke = ceil((ionBeamWidth + subLength) / ionBeamWidth);
upRange = nUpStroke * ionBeamWidth;
nDownStroke = ceil((ionBeamWidth/ms + subLength) / ionBeamWidth);
downRange = nDownStroke * ionBeamWidth;

preSize = [0, ionBeamWidth-(ionBeamWidth/ms)];
zPadPre = 0.5 * padarray(z, preSize, nan, 'pre');
postSize = [0, upRange-(ionBeamWidth+subLength)];
zPadPost = padarray(zPadPre, postSize, nan, 'post');

zUp = padarray(zPadPost, [leafWidth/2, 0], nan, 'both');

postSize2 = [0, downRange-(ionBeamWidth/ms + subLength)];
zPadPost2 = 0.5 * padarray(z, postSize2, nan, 'post');

zDown = padarray(zPadPost2, [leafWidth/2, 0], nan, 'both');

% figure, mesh(zUp);
% axis equal;
% figure, mesh(zDown);
% axis equal;

xform = [ 1    0   0
         -shx  1   0
          0    0   1 ];
tform = maketform('affine', xform);
zr = imtransform(zUp, tform, 'FillValues', nan);  % rising strokes
% figure, mesh(zr);
% axis equal;

xform2 = [ 1   0   0
          shx  1   0
           0   0   1 ];
tform2 = maketform('affine',xform2);
zf = imtransform(zDown, tform2, 'FillValues', nan);  % falling strokes
% figure, mesh(zf);
% axis equal;

% divide rising strokes
startIx = 1:ionBeamWidth:upRange-ionBeamWidth+1;
endIx = [startIx(2:end)-1, upRange];

upRibbons = nan(size(zr,1), nUpStroke);

for i = 1:nUpStroke
    upRibbons(:, i) = mean(zr(:, startIx(i):endIx(i)), 2);
end

figure;
hUp = ribbon(upRibbons);
xlabel('X')
ylabel('Y')
zlabel('Z')

% figure, waterfall(upRibbons');

% divide falling strokes
startIx = 1:ionBeamWidth:downRange-ionBeamWidth+1;
endIx = [startIx(2:end)-1, downRange];

downRibbons = nan(size(zf,1), nDownStroke);

for i = 1:nDownStroke
    downRibbons(:, i) = mean(zf(:, startIx(i):endIx(i)), 2);
end

figure;
hDown = ribbon(downRibbons);
xlabel('X')
ylabel('Y')
zlabel('Z')

% figure, waterfall(downRibbons');

ribbons = [upRibbons downRibbons];
figure, hold on;
plot(ribbons)

% meanEtchTime = nan(vStroke, 1);
% for j = 1:vStroke
%     row_j = ribbons(j, :);
%     row_j = row_j(~isnan(row_j));
%     num_j = numel(row_j);
%     if num_j > 0
%         meanEtchTime(j) = sum(row_j) / num_j;
%     end
% end
% 
% plot(meanEtchTime, 'k--')

ribbonNos = ribbons(~isnan(ribbons));
meanEtchTimeG = mean(ribbonNos); % global mean etch time

plot(meanEtchTimeG*ones(vStroke, 1), 'k--')

xlim([1 vStroke])

maxEtchTimeG = max(ribbonNos);
minEtchTimeG = min(ribbonNos);

% max tunable ratio of dwell time
maxTunableRatio = leafWidth/vStroke;

% maxEtchTime < strokeTime < minEtchTime/(1-maxTunableRatio)
if maxEtchTimeG < minEtchTimeG/(1-maxTunableRatio)
    strokeTimeRange = [maxEtchTimeG, minEtchTimeG/(1-maxTunableRatio)];
    strokeTimeOpt = meanEtchTimeG / (1 - 0.5*maxTunableRatio);
    if strokeTimeOpt>strokeTimeRange(1) && strokeTimeOpt<strokeTimeRange(2)
        strokeTime = strokeTimeOpt;
    end
end

plot(maxEtchTimeG*ones(vStroke, 1), 'k:');
plot(minEtchTimeG/(1-maxTunableRatio)*ones(vStroke, 1), 'k:');

plot(strokeTime*ones(vStroke, 1), 'r-')
plot((2*meanEtchTimeG-strokeTime)*ones(vStroke, 1), 'b:')

elapsedTime = ms*strokeTime;

%
upDwells = strokeTime - upRibbons;
downDwells = strokeTime - downRibbons;

figure;
plot(upDwells)
