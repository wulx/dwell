%% tansform and slice dwell-time map
close all; clear all; clc

% load etch depth data
load etch_depth2.mat
z = depth; % in nm
x = 0.5:399.5; % mm, 400 points in total
y = 0.5:399.5; % mm
figure, meshc(x, y, z);
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
zr = imtransform(zUp, tform, 'nearest', 'FillValues', nan);  % rising strokes
% figure, mesh(zr);
% axis equal;

xform2 = [ 1   0   0
          shx  1   0
           0   0   1 ];
tform2 = maketform('affine',xform2);
zf = imtransform(zDown, tform2, 'nearest', 'FillValues', nan);  % falling strokes
% figure, mesh(zf);
% axis equal;

% divide rising strokes
rStartIx = 1:ionBeamWidth:upRange-ionBeamWidth+1;
rEndIx = [rStartIx(2:end)-1, upRange];

upRibbons = nan(size(zr,1), nUpStroke);

for i = 1:nUpStroke
    zri = zr(:, rStartIx(i):rEndIx(i));
    for ii = 1:size(zri,1)
        zrii = zri(ii, :);
        idx = ~isnan(zrii);
        c = sum(idx);
        if c > 0
            upRibbons(ii, i) = sum(zrii(idx)) / c;
        end
    end
end

% figure;
% hUp = ribbon(upRibbons);
% xlabel('X')
% ylabel('Y')
% zlabel('Z')

zr2 = zr;
zr2(:, rStartIx) = nan;
zr2 = imtransform(zr2, tform2, 'FillValues', nan);
figure, mesh(zr2);
axis equal

figure, waterfall(upRibbons');

% divide falling strokes
fStartIx = 1:ionBeamWidth:downRange-ionBeamWidth+1;
fEndIx = [fStartIx(2:end)-1, downRange];

downRibbons = nan(size(zf,1), nDownStroke);

for i = 1:nDownStroke
    zfi = zf(:, fStartIx(i):fEndIx(i));
    for ii = 1:size(zri,1)
        zfii = zfi(ii, :);
        idx = ~isnan(zfii);
        c = sum(idx);
        if c > 0
            downRibbons(ii, i) = sum(zfii(idx)) / c;
        end
    end
end

% figure;
% hDown = ribbon(downRibbons);
% xlabel('X')
% ylabel('Y')
% zlabel('Z')

zf2 = zf;
zf2(:, fStartIx) = nan;
zf2 = imtransform(zf2, tform, 'FillValues', nan);
figure, mesh(zf2);
axis equal

figure, waterfall(downRibbons');

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
upDwells = upDwells(leafWidth/2+1:end-leafWidth/2, :);
downDwells = strokeTime - downRibbons;
downDwells = downDwells(leafWidth/2+1:end-leafWidth/2, :);

figure, plot(upDwells);
figure, plot(downDwells);

