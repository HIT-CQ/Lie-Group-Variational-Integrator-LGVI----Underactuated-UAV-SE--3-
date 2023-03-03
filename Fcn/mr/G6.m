function G_rho = G6(expc6)
% d老 = G(老)*缶  
% Input:    硌杅釴梓 老=[成;汕]  (6℅1)
% Output: 撻淝 G(老)
%    
Theta = expc6(1:3);     % 成
y = norm(Theta);          % 牟

A1 = 2/(y^2) - 3/(4*y)*cot(y/2) - 1/8*(csc(y/2))^2;
A2 = 1/(y^4) - 1/(4*y^3)*cot(y/2) - 1/(8*y^2)*(csc(y/2))^2;

G_rho = eye(6) + 0.5*ad(expc6) + A1*(ad(expc6))^2 + A2*(ad(expc6))^4;    
end