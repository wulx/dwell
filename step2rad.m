function radians = step2rad(nSteps, stepAngleDeg)
%STEP2RAD convert number of steps of the specified stepper motor to radians

if any(~isfinite(nSteps) | nSteps<0)
    error('stepper numbers should be finite non-negative')
end

if ~isfinite(stepAngleDeg) || stepAngleDeg<0
    error('step angle should be finite non-negative')
end

radians = deg2rad(stepAngleDeg * nSteps);

end