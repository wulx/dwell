%PROFOPT motion profile optimization to obtain the minimum RMSE of dwell time
clear all; close all; clc

%% load objective data
load dwell_time.mat

for n = 1:8
    ogee = upDwellTime(:, n);
    
    [xx, ys] = spis(ogee);
    
    % low pass filtering
    lowPass = ys>maxDwellTime;
    if any(lowPass)
        ys(lowPass) = maxDwellTime;
    end
    
    figure, hold on;
    plot(xx, ys, 'k-')
    
    opts = {'MINPEAKHEIGHT', -1, 'MINPEAKDISTANCE', 15, 'NPEAKS', 5};
    [locs, ptype] = findpoints(ys, opts);
    
    % critical points
    for i = 1:numel(locs)
        li = locs(i);
        ti = ptype(i);
        if ti > 0 % peaks
            plot(xx(li), ys(li), 'ro')
        elseif ti < 0 % valleys
            plot(xx(li), ys(li), 'bo')
        else % inflection points
            plot(xx(li), ys(li), 'ko')
        end
    end
    
    % boundaries
    plot(xx([1 end]), ys([1 end]), 'go')
    
    % set params
    Locs = [1 locs numel(xx)];
    dwells = ys(Locs);
    angles = asin(dwells/maxDwellTime);
    degs = rad2deg(angles);
    
    
    %ionWidth = 60; % ion beam width in mm
    %vTrans = ionWidth / elapsedTime;
    
    stroke = 460;
    
    stepAngleDeg = 1.8 / 8;
    

    stps = round(degs / stepAngleDeg);
    
    dstps = abs(stps(2:end) - stps(1:end-1));
    
    lowStp = dstps<4;
    lowStp2 = [lowStp(1) lowStp(1:end)] | [lowStp(1:end) lowStp(end)];
    
    amp = 3*lowStp2;
    errbar = [0 ptype 0];

    stps = stps + amp.*errbar + 3*errbar;
    
    
    S = cell(1, numel(Locs)-1);
    T = cell(size(S));
    for j = 2:numel(degs)
        
        sn_tot = (stps(j) - stps(j-1));
        
        
        strokeDown = 0;
        if sn_tot < 0
            sn_tot = -sn_tot;
            strokeDown = 1;
            disp('falling down')
        end
        
        if sn_tot > 60
            s_u = 2;
        else
            s_u = 1;
        end
        
        
        sn_tot = sn_tot/s_u;
        t_tot = strokeTime * (xx(Locs(j)) - xx(Locs(j-1))) / stroke;
        

        w_a = 1/3; %<--------------------------------------------------%
        sn_a = round(w_a * sn_tot);
        w_d = 1/3;  %<-----------------------------------------------------------%
        %     if w_d > (1-w_a), w_d = (1-w_a)*w_d; end
        
        sn_d = round(w_d * sn_tot);
        sn_c = sn_tot - sn_a - sn_d;
        
        a_max = 60000; % benchtest result: 66333 pulse/s^2
        
        % 0 <= w_f <= 1
        w_f = rand * (sn_tot/t_tot - 2)/(0.25*a_max*t_tot);  %<-------------------%
        
        % (0, (sn_tot/t_tot - 2)/(0.25*a_max*t_tot)]
        f_i = max(floor(sn_tot/t_tot - 0.25*a_max*t_tot*w_f), 2);
        
        f_m = sym('f_m', 'positive');
        
        % solve multiple equations of linear speed ramp
        sols = solve(2*sn_a/(f_i+f_m) + sn_c/f_m + 2*sn_d/(f_i+f_m) == t_tot, ...
            f_m^2 - f_i^2 <= 2*a_max*min(sn_a, sn_d), ...
            f_m >= sn_tot/t_tot, f_m);
        
        f_m = round( double(sols) );
        
%         if ~isempty(f_m)
%             2*(sn_a+sn_d)/(f_i+f_m) + sn_c/f_m
%         end
        
        % range of initial frequency
        
        
        
        sn = [sn_a, sn_c, sn_d];
        pf = round([f_i, f_m]);
        method = 'round';
        [f_list, dt_list] = time_per_step(sn, pf, s_u, method);
        
        %sn_plot(f_list, dt_list)
        %title(num2str( [sum(dt_list) sn pf]))
        
        % #2 step time -----------------------------------------------------------%
        timeSeqs = steptime(f_list, dt_list);
        
        steps = 0:numel(timeSeqs);
        
        timeline = arrayfun(@(t) sum(timeSeqs(1:t)), steps);
        
        if strokeDown
            steps = stps(j-1) - steps;
        else
            steps = stps(j-1) + steps;
        end
%         figure, stairs(timeline, steps);
        
        T{j-1} = timeSeqs;
        if j < numel(degs)
            S{j-1} = steps(1:end-1);
        else
            S{j-1} = steps;
        end
        
    end
    
    % cat
    T = cell2mat(T);
    timeline = [0 arrayfun(@(t) sum(T(1:t)), 1:numel(T))];
    S = cell2mat(S);
%     stairs(timeline, S);
    
    timeStep = 0.001;
    totalTime = sum( T(:) );
    nScan = ceil(totalTime / timeStep);
    
    steps = zeros(1, nScan);
    timeline = linspace(0, totalTime, nScan);
    
    nStep = numel(T);
    for i = 1:nStep
        steps(timeline>=sum(T(1:i-1)) & timeline<=sum(T(1:i))) = S(i);
    end
    
%     figure, plot(steps)
%     
    
    leafWidth = 60;
    initAngleDeg = 0;
    projWidths = step2width(steps, stepAngleDeg, leafWidth, initAngleDeg);
    
    scaleDivs = xx;
    dwellTime = timecount(projWidths, strokeTime, scaleDivs);
%     figure, plot(scaleDivs, dwellTime);
    
    filt = ~[1 isnan(ogee)' 1];
    scaleDivs = scaleDivs(filt);
    dwellTime = dwellTime(filt);
    
    figure, hold on;
    plot(scaleDivs, dwellTime, 'r-');
    isNum = ~isnan(ogee);
    plot(scaleDivs, ogee(isNum), 'k-')
    
end


%% outputs
% root-mean-square deviation of dwell time
% r = rmse(dwellTime, dwellTimeEstimate);