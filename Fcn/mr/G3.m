function G_Theta = G3(expc3) 
% Input:  姿态指数坐标Θ (3×1向量)         
% Output: G(Θ)  
%   dΘ = G(Θ)*omega

theta = norm(expc3);    % θ = ||Θ||

G_Theta = eye(3) + 1/2*(VecToso3(expc3)) + ...
         (1/((theta)^2) - (1+cos(theta))/(2*theta*sin(theta)))*(VecToso3(expc3))^2;
    
end