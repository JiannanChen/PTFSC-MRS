function [sys,x0,str,ts,simStateCompliance] = mrs_ptfsc(t,x,u,flag)

%   Copyright 1990-2010 The MathWorks, Inc.

%
% The following outlines the general structure of an S-function.
%
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=[];

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=[];

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=[];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 1;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 9;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = 5.5;

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,~)
% MRS物理参数 - MRS physical paras
M = 200;  
Jw = 0.005; 
c = 0.05; 
r = 0.1;
k = 6; 

% 系统参数  - system paras
a1 = 2*c/(M*r^2 + 2*Jw); 
b1 = k*r/(M*r^2 + 2*Jw);

% 控制参数
delta = 0.1; 
epsilon = 0.1; 
eta = 0.6;
m = 1 - 0.5*eta;
n = 1 + 0.5*eta;

% T = 2; Td = 0.5;
T = 1; Td = 0.1;
k1 = (pi/(eta*T))*0.5*3^(eta/2);
% k1 = 100;

% 期望跟踪信号 - desired signals
x1d = 3*sin(t); 
dx1d = 3*cos(t); 

% 误差信号 
z1 = x(1) - x1d; 

if t >= Td
    rho = 1;
else
    rho = 1 - (1-t/Td)*exp(1-(Td/(Td-t)));
end
    
if t >= Td
    drho = 0;
else
    drho = (1/Td) *exp(1-(Td/(Td-t)))  - (1-t/Td)*exp(1-(Td/(Td-t))) * (-Td/(Td-t)^2) ;
end
      
z1hat = rho * z1; 

% 限制函数 - constraint function
kb1 = 5 + 0.1*sin(t); 
dkb1 = 0.1*cos(t); 

% 控制方案 - control scheme
R1 = 1 / (kb1^2 - z1hat^2); 

if z1hat^4 > delta^4
    Phi1 = 1;
else
    Phi1 = sin(0.5*pi*(z1hat^4 / delta^4)) * sin(0.5*pi*(z1hat^4 / delta^4));
end

% 控制输入 - control input 
if z1hat^4 > delta^4  % 针对Phi1*z1hat^(2*m-2)的近似，sinx 约等于 x
    var1 = sign(z1hat)*abs(z1hat)^(2*m-2);
else
    var1 = (0.5*pi*(z1hat^4 / delta^4)) * (0.5*pi*(z1hat^2 / delta^4))*sign(z1hat)*(abs(z1hat)^(2*m));
end

% uv = (1/b1)*(a1*x(1) + dx1d - epsilon*R1*drho^2*z1^3 - epsilon*R1*(dkb1/kb1)^2*z1*z1hat^2 - k1*z1...
%     - k1*Phi1*z1hat^(2*m-2)*R1^(m-1)*z1 - k1*Phi1*R1^(n-1)*z1hat^(2*n-2)*z1);
uv = (1/b1)*( a1*x(1) + dx1d - epsilon*R1*(drho^2)*z1^3 - epsilon*R1*((dkb1/kb1)^2)*z1*z1hat^2 - k1*z1...
    - k1*var1*R1^(m-1)*z1 - k1*Phi1*(R1^(n-1))*sign(z1hat)*(abs(z1hat)^(2*n-2))*z1 );

% 系统动态模型 - system dynamic model
sys(1) = -a1*x(1) + b1*uv;

% end mdlDerivatives

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,~)
% MRS物理参数 - MRS physical paras
M = 200;  
Jw = 0.005; 
c = 0.05; 
r = 0.1;
k = 6; 

% 系统参数  - system paras
a1 = 2*c/(M*r^2 + 2*Jw); 
b1 = k*r/(M*r^2 + 2*Jw);

% 控制参数
delta = 0.1; 
epsilon = 0.1; 
eta = 0.6;
m = 1 - 0.5*eta;
n = 1 + 0.5*eta;

% T = 2; Td = 0.5;
T = 1; Td = 0.1;
k1 = (pi/(eta*T))*0.5*3^(eta/2);
% k1 = 100;

% 期望跟踪信号 - desired signals
x1d = 3*sin(t); 
dx1d = 3*cos(t); 

% 误差信号 
z1 = x(1) - x1d; 

if t >= Td
    rho = 1;
else
    rho = 1 - (1-t/Td)*exp(1-(Td/(Td-t)));
end
    
if t >= Td
    drho = 0;
else
    drho = (1/Td) *exp(1-(Td/(Td-t)))  - (1-t/Td)*exp(1-(Td/(Td-t))) * (-Td/(Td-t)^2) ;
end
      
z1hat = rho * z1; 

% 限制函数 - constraint function
kb1 = 5 + 0.1*sin(t); 
dkb1 = 0.1*cos(t); 

% 控制方案 - control scheme
R1 = 1 / (kb1^2 - z1hat^2); 

if z1hat^4 > delta^4
    Phi1 = 1;
else
    Phi1 = sin(0.5*pi*(z1hat^4 / delta^4)) * sin(0.5*pi*(z1hat^4 / delta^4));
end

% 控制输入 - control input 
if z1hat^4 > delta^4  % 针对Phi1*z1hat^(2*m-2)的近似，sinx 约等于 x
    var1 = sign(z1hat)*abs(z1hat)^(2*m-2);
else
    var1 = (0.5*pi*(z1hat^4 / delta^4)) * (0.5*pi*(z1hat^2 / delta^4))*sign(z1hat)*(abs(z1hat)^(2*m));
end

% uv = (1/b1)*(a1*x(1) + dx1d - epsilon*R1*drho^2*z1^3 - epsilon*R1*(dkb1/kb1)^2*z1*z1hat^2 - k1*z1...
%     - k1*Phi1*z1hat^(2*m-2)*R1^(m-1)*z1 - k1*Phi1*R1^(n-1)*z1hat^(2*n-2)*z1);
uv = (1/b1)*( a1*x(1) + dx1d - epsilon*R1*(drho^2)*z1^3 - epsilon*R1*((dkb1/kb1)^2)*z1*z1hat^2 - k1*z1...
    - k1*var1*R1^(m-1)*z1 - k1*Phi1*(R1^(n-1))*sign(z1hat)*(abs(z1hat)^(2*n-2))*z1 );

% 系统状态 - system state
sys(1) = x(1);
sys(2) = z1;
sys(3) = z1hat;
% 其余信息 - other info
sys(4) = Phi1;
sys(5) = R1;
sys(6) = rho;
sys(7) = x1d;
sys(8) = kb1;
sys(9) = -kb1;

% end mdlOutputs
