function [bd,vd,dvd] = desired_trans(t)

% Input:     t
% Output:  当前时刻t时的期望轨迹

bd = [0.4*sin(pi*t);0.6*cos(pi*t);0.4*t];
vd = [0.4*pi*cos(pi*t);-0.6*pi*sin(pi*t);0.4];
dvd = [-0.4*pi^2*sin(pi*t);-0.6*pi^2*cos(pi*t);0];

end
