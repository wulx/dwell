% function setramp
close all; clear all; clc

s_u = 2;
t_tot = 4 / s_u;
sn_tot = 168;
w_a = max(rand, 0.1); %<--------------------------------------------------%
sn_a = round(w_a * sn_tot);
w_d = rand;  %<-----------------------------------------------------------%
if w_d > (1-w_a), w_d = (1-w_a)*w_d; end

sn_d = round(w_d * sn_tot);
sn_c = sn_tot - sn_a - sn_d;

a_max = 10000;

% 0 <= w_f <= 1
w_f = rand * (sn_tot/t_tot - 2)/(0.25*a_max*t_tot);  %<-------------------%

% (0, (sn_tot/t_tot - 2)/(0.25*a_max*t_tot)]
f_i = round( sn_tot/t_tot - 0.25*a_max*t_tot*w_f );

f_m = sym('f_m', 'positive');

% solve multiple equations of linear speed ramp
sols = solve(2*sn_a/(f_i+f_m) + sn_c/f_m + 2*sn_d/(f_i+f_m) == t_tot, ...
    f_m^2 - f_i^2 <= 2*a_max*min(sn_a, sn_d), ...
    f_m >= sn_tot/t_tot, f_m);

f_m = round( double(sols) );

if ~isempty(f_m)
    2*(sn_a+sn_d)/(f_i+f_m) + sn_c/f_m
end

% range of initial frequency

sn = [sn_a, sn_c, sn_d];
pf = round([f_i, f_m]);
method = 'round';
s_u = 1;
[f_list, dt_list] = time_per_step(sn, pf, s_u, method);

sn_plot(f_list, dt_list)
title(num2str( [sum(dt_list) sn pf]))
