function [f]=trans_control(R,bt,nu,vd,dvd)

global m g e3 P L

f = e3'*R'*(m*g*e3+P*bt+L*(R*nu-vd)-m*dvd);

end
