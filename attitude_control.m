function [tau]=attitude_control(Om,Omd,dOmd,Q)

global J e1 e2 e3 La p ki kp a1 a2 a3

om=Om-Q.'*Omd;

s = a1*cross(Q.'*e1,e1)+a2*cross(Q.'*e2,e2)+a3*cross(Q.'*e3,e3);
H = eye(3,3)-((2*(1-(1/p)))/(s.'*s))*(s*s.');
w = a1*cross(e1,cross(om,Q.'*e1))+a2*cross(e2,cross(om,Q.'*e2))+a3*cross(e3,cross(om,Q.'*e3));
z = s/((s.'*s)^(1-(1/p)));
Psi = om+ki*z;

%dOmd=dOm+(ki/((s'*s)^(1-(1/p))))*H*w;

tau =J*(Q.'*dOmd-(ki*H/(s.'*s)^(1-(1/p)))*w)...
     +VecToso3(Q.'*Omd)*J*(Q.'*Omd-ki*z)+ki*J*(cross(z,(Q.'*Omd)))...
     +cross(ki*J*(w+Q.'*Omd),z)-kp*s...
     -(La*Psi/(Psi.'*La*Psi)^(1-(1/p)));


