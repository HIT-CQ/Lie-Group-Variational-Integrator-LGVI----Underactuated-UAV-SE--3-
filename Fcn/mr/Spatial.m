function Spatial_Inertia_Matrix = Spatial(Inertia,mass)
% Input:    ���Ծ���3��3    &    ����
% Output: �ռ��������(Spatial Inertia Matrix)  6��6
   
Spatial_Inertia_Matrix = [Inertia        zeros(3);
                                        zeros(3)     mass*eye(3)];
end

