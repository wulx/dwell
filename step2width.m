function projWidths = step2width(steps, stepAngleDeg, leafWidth, initAngleDeg)
%STEPWIDTH convert number of steps of the specified stepper motor to
%projected widths

if any(~isfinite(steps) | steps<0)
    error('number of steps should be finite non-negative')
end

if ~isfinite(stepAngleDeg) || stepAngleDeg<0
    error('step angle should be finite non-negative')
end

rotAngles = deg2rad(initAngleDeg + stepAngleDeg*steps);

% the scanning leaf parallels to ion beam when rotation angle is 0,
% and perpendicular when the angle is  pi/2.
projWidths = abs(cos(rotAngles) * leafWidth);

end