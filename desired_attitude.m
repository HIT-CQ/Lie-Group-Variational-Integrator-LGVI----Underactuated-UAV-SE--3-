function [Rd] = desired_attitude(bt,vt,dvd)

% Input:   bt,vt,dvd  
% Output:  根据位置轨迹 输出 期望姿态

global m g e1 e3 P L mu

r3d = (m*g*e3+P*bt+L*vt-m*dvd)/norm(m*g*e3+P*bt+L*vt-m*dvd);     % 式(19)
sd   = cross(ones(3,1),r3d)+mu*e1;
r2d = cross(r3d,sd)/norm(cross(r3d,sd));
Rd  =[cross(r2d,r3d) r2d r3d];

end
