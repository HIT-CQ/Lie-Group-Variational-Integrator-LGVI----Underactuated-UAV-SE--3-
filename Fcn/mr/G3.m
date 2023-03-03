function G_Theta = G3(expc3) 
% Input:  ��ָ̬�����ꦨ (3��1����)         
% Output: G(��)  
%   d�� = G(��)*omega

theta = norm(expc3);    % �� = ||��||

G_Theta = eye(3) + 1/2*(VecToso3(expc3)) + ...
         (1/((theta)^2) - (1+cos(theta))/(2*theta*sin(theta)))*(VecToso3(expc3))^2;
    
end