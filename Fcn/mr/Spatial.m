function Spatial_Inertia_Matrix = Spatial(Inertia,mass)
% Input:    惯性矩阵3×3    &    质量
% Output: 空间惯量矩阵(Spatial Inertia Matrix)  6×6
   
Spatial_Inertia_Matrix = [Inertia        zeros(3);
                                        zeros(3)     mass*eye(3)];
end

