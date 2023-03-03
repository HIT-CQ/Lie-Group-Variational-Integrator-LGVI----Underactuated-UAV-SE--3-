function AdTInv = AdjointInv(T)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes T a transformation matrix SE3. 
% Returns the corresponding 6x6 adjoint representation Inverse [AdT]^-1.
% [AdT]^-1 = [AdT^-1] = [RT  0; -RT*VecToso3(p)  RT]
% Example Input:
% 
% clear; clc;
% T = [[1, 0, 0, 0]; [0, 0, -1, 0]; [0, 1, 0, 3]; [0, 0, 0, 1]];
% AdT = Adjoint(T)
% 
% Output:
% AdT =
%      1     0     0     0     0     0
%      0     0     1     0     0     0
%      0    -1     0     0     0     0
%      0     3     0     1     0     0
%      0     0     0     0     0     1
%      3     0     0     0    -1     0

[R, p] = TransToRp(T);
AdTInv = [transpose(R), zeros(3); -transpose(R)*VecToso3(p), transpose(R)];
end