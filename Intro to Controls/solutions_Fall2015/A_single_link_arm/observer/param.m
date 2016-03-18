clear all

%--------------------------------------------------------
% actual parameters - used in physics model in simulation
% actual system parameters
AP.m = 0.5;  % kg
AP.ell = 0.3; % m
AP.b = 0.01; % N m s
AP.g = 9.8; % m/s^2
% actual initial conditions
AP.theta0 = 0;
AP.thetadot0 = 0;
% input constraint
AP.tau_max = 1;

%--------------------------------------------------------
% control parameters - these are the parameters known to the controller,
% which may not be the actual parameters
P.m = AP.m;%*(0.9);  % kg
P.ell = AP.ell;%*(0.95); % m
P.b = AP.b;%*(1.1); % N m s
P.g = AP.g; % m/s^2
% initial conditions
P.theta0 = AP.theta0;
P.thetadot0 = AP.thetadot0;
% input constraint
P.tau_max = AP.tau_max;

% equalibrium torque
P.theta_e = 0;
P.tau_e = P.m*P.g*P.ell/2*cos(P.theta_e);
P.tau_sat = P.tau_max-P.tau_e;

% select PD gains
A_th = 50*pi/180;
zeta = 0.707;
P.kp = (P.tau_sat)/A_th;
wn = sqrt(3*P.kp/P.m/P.ell^2);
P.kd = P.m*P.ell^2/3*2*zeta*wn-P.b;

kp = (P.tau_max-P.tau_e)/A_th;
wn = sqrt(3*kp/P.m/P.ell^2);
kd = P.m*P.ell^2/3*2*zeta*wn-P.b;

% sample rate
P.Ts = 0.01;

% state space model
P.A = [0 1; 0 -3*P.b/P.m/(P.ell)^2];
P.B = [0;3/P.m/(P.ell^2)];
P.C = [1 0];
P.D = 0;

% check the state space model to make sure it is correct
%[A1,B1,C1,D1] = linmod('arm_sim_lin');

% augmented state space model to include integrator
AA = [P.A, zeros(2,1); P.C, 0];
BB = [P.B; 0];

% desired closed loop pole locations
p_cl = [roots([1,2*zeta*wn,wn^2]);-.1];
tmp = place(AA,BB,p_cl);
P.K = tmp(:,1:2);
P.kr = -1/(P.C*inv(P.A-P.B*P.K)*P.B);
P.ki = tmp(1,3);

% observer poles
p_obs = 10*p_cl(1:end-1);
P.L = place(P.A',P.C',p_obs)';
eig(P.A-P.B*P.K-P.L*P.C)


% ekf parameters
P.Q = diag([.000000001,1]);
P.R = .0001;


