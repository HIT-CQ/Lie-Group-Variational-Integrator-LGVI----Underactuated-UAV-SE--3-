%% 文献：
% Viswanathan S P, Sanyal A K, Samiei E. 
% Integrated guidance and feedback control of underactuated robotics system in SE (3)
% Journal of Intelligent & Robotic Systems
% 2018, 89: 251-263

% 代码参考 https://github.com/ssingh130/LGVI_FTS

%% 加载函数路径
addpath(genpath(pwd));

%% 参数
global J M m e1 e2 e3 g mu...
    P L La p ki kp a1 a2 a3

% UAV参数
J = diag([0.0820,0.0845,0.1377]);
m = 4.34;                                 
M = m*eye(3);                             
g = 9.81;
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1]; 
mu = 4;

% 控制增益参数
P = 38*eye(3);
L = 25*eye(3);
La = 3.5*eye(3);
p = 0.75;
ki = 0.04;
kp = 4.5;
a1 = 1.2;
a2 = 1.1;
a3 = 1;
% A = diag([a1,a2,a3]);

% 初始状态参数
b0 = [1;0;0];        % 初始位矢
R0 = eye(3);        % 初始姿态
Om0=[0;0;0];      % 初始角速度
nu0 = [0;0;0];      % 初始线速度(体坐标系)

%% 仿真时间
h = .01;     % 步长
t0 = 0;
tf = 5;
TimeLength = tf - t0;
n = TimeLength/h;
t = linspace(t0,tf,n);


%% 预分配内存
Rd       = zeros(3,3,length(t));
bd       = zeros(3,length(t));
vd       = zeros(3,length(t));
dvd     = zeros(3,length(t));
Omd   = zeros(3,length(t));
dOmd = zeros(3,length(t));

R         = zeros(3,3,length(t));
b         = zeros(3,length(t));
Om     = zeros(3,length(t));
nu       = zeros(3,length(t));
v         = zeros(3,length(t));

bt        = zeros(3,length(t));
vt        = zeros(3,length(t));
Q        = zeros(3,3,length(t));

F         = zeros(3,3,length(t));
Fd       = zeros(3,3,length(t));

f          = zeros(1,length(t));
tau      = zeros(3,length(t));

%% 初始赋值
[bd(:,1),vd(:,1),dvd(:,1)] = desired_trans(t(1));

R(:,:,1)    = R0;
Om(:,1) = Om0;
Omd(:,1) = Om0;
b(:,1)     = b0;
bt(:,1)    = b(:,1) - bd(:,1);
nu(:,1)    = nu0;
v(:,1)=  R(:,:,1)*nu(:,1);
vt(:,1)     = R(:,:,1)*nu(:,1) - vd(:,1);
Rd(:,:,1)  = desired_attitude(bt(:,1),vt(:,1),dvd(:,1));

%% LGVI (李群变分积分器)
for k = 1:n-1
    
    [bd(:,k+1),vd(:,k+1),dvd(:,k+1)] = desired_trans(t(k+1));
    
    F(:,:,k)      = MatrixExp3(h*VecToso3(Om(:,k)));
    R(:,:,k+1) = R(:,:,k) * F(:,:,k);
    b(:,k+1)   = h*R(:,:,k)*nu(:,k)+b(:,k);
    
    % 控制推力
    bt(:,k+1) = b(:,k+1)-bd(:,k+1);
    f(:,k)   = trans_control(R(:,:,k),bt(:,k),nu(:,k),vd(:,k),dvd(:,k));
    %
    
    nu(:,k+1) = M\( m*F(:,:,k).'*nu(:,k) + h*m*g*R(:,:,k+1).'*e3 - h*f(:,k)*e3 );
    v(:,k+1)=R(:,:,k+1)*nu(:,k+1);
    vt(:,k+1)=v(:,k+1)-vd(:,k+1);
    
    Rd(:,:,k+1) = desired_attitude(bt(:,k+1),vt(:,k+1),dvd(:,k+1));
    Fd(:,:,k+1)=Rd(:,:,k).'*Rd(:,:,k+1);
    Omd(:,k+1) = so3ToVec(1/h*MatrixLog3(Fd(:,:,k+1)));
    dOmd(:,k)= (Omd(:,k+1)-Omd(:,k))/h;
    Q(:,:,k)=Rd(:,:,k)'*R(:,:,k);
    
    % 控制力矩
    tau(:,k) = attitude_control(Om(:,k),Omd(:,k),dOmd(:,k),Q(:,:,k));
    %
    
    Om(:,k+1) = J\( F(:,:,k)'*J*Om(:,k)+h*tau(:,k) );
    
end

Q(:,:,n)=Rd(:,:,n)'*R(:,:,n);
f(:,n)   = trans_control(R(:,:,n),bt(:,n),nu(:,n),vd(:,k),dvd(:,n));
tau(:,n) =attitude_control(Om(:,n),Omd(:,n),dOmd(:,n),Q(:,:,n));

%% Plot
figure('name','trans_error')
plot(t,bt);