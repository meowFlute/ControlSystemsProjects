function [sys,x0,str,ts,simStateCompliance] = mav_dynamics(t,x,u,flag,P)

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(P);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u,P);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);

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
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(P)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 12;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 12;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [...
    P.pn0;...
    P.pe0;...
    P.pd0;...
    P.u0;...
    P.v0;...
    P.w0;...
    P.phi0;...
    P.theta0;...
    P.psi0;...
    P.p0;...
    P.q0;...
    P.r0;...
    ];

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
%%
function sys=mdlDerivatives(t,x,uu, P)

    pn    = x(1);
    pe    = x(2);
    pe    = x(3);
    u     = x(4);
    v     = x(5);
    w     = x(6);
    phi   = x(7);
    theta = x(8);
    psi   = x(9);
    p     = x(10);
    q     = x(11);
    r     = x(12);
    fx    = uu(1);
    fy    = uu(2);
    fz    = uu(3);
    ell   = uu(4);
    m     = uu(5);
    n     = uu(6);
    
    %% build the rotation matrix
    rot1 = [cos(psi)    sin(psi)    0;...
           -sin(psi)    cos(psi)    0;...
            0           0           1;...
            ];
    rot2 = [cos(theta)  0          -sin(theta);...
            0           1           0;...
            sin(theta)  0           cos(theta);...
            ];
    rot3 = [1           0           0;...
            0           cos(phi)    sin(phi);...
            0          -sin(phi)    cos(phi);...
            ];
    Rv_b = rot3*rot2*rot1; % !!! From the vehicle to the body frame !!!
    
    %%
    position_dot = Rv_b' * [u; v; w;];       % See eqn 3.14 or 3.1
    pndot       = position_dot(1);      % North position differential equation
    pedot       = position_dot(2);      % East position differential equation
    pddot       = position_dot(3);      % Down position differential equation
    
    udot        = (r*v) - (q*w) + ((1/P.mass)*fx);  % X-axis velocity differential equation
    vdot        = (p*w) - (r*u) + ((1/P.mass)*fy);  % Y-axis velocity differential equation
    wdot        = (q*u) - (p*v) + ((1/P.mass)*fz);  % Z-axis velocity differential equation
    
    %% Rotational Kinematics
    R = [1      sin(phi)*tan(theta)     cos(phi)*tan(theta);...
        0           cos(phi)                -sin(phi);...
        0       sin(phi)*sec(theta)     cos(phi)*sec(theta)];
    rotation_dot = R*[p; q; r];
    %differential equations 
    phidot = rotation_dot(1);
    thetadot = rotation_dot(2);
    psidot = rotation_dot(3);
    
    Gamma = P.Jx*P.Jz - P.Jxz^2;
    Gamma_1 = (P.Jxz*(P.Jx- P.Jy + P.Jz))/Gamma;
    Gamma_2 = (P.Jz*(P.Jz-P.Jy) + P.Jxz^2)/Gamma;
    Gamma_3 = P.Jz/Gamma;
    Gamma_4 = P.Jxz/Gamma;
    Gamma_5 = (P.Jz-P.Jx)/P.Jy;
    Gamma_6 = P.Jxz/P.Jy;
    Gamma_7 = ((P.Jx-P.Jy)*P.Jx + P.Jxz^2)/Gamma;
    Gamma_8 = P.Jx/Gamma;
    
    pdot = Gamma_1*p*q - Gamma_2*q*r + Gamma_3*ell + Gamma_4*n;
    qdot = Gamma_5*p*r - Gamma_6*(p^2 - r^2) + (1/P.Jy)*m;
    rdot = Gamma_7*p*q - Gamma_1*q*r + Gamma_4*ell + Gamma_8*n;

sys = [pndot; pedot; pddot; udot; vdot; wdot; phidot; thetadot; psidot; pdot; qdot; rdot];

% end mdlDerivatives
%%
%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

sys = x;

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
