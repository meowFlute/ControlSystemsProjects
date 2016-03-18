clear;
% gain on dirty derivative
P.tau = 5;
P.gravity = 9.8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params for sensors from appendix H
P.gyro_max          = 350*pi/180; % + or - 350 degrees/s
gyro_bandwidth      = 80; %Hz
gyro_noise_density  = 0.015*pi/180; % rad/s/sqrt(Hz)
P.sigma_gyro        = sqrt(gyro_bandwidth)*gyro_noise_density; % standard deviation of error in deg/s

P.accel_max         = 6*P.gravity; % + or - 6g's
accel_bandwidth     = 100; % Hz
accel_noise_density = 250e-6*P.gravity; %250 micro gravities/sqrt(Hz)
P.sigma_accel       = accel_noise_density*sqrt(accel_bandwidth);

P.beta_static       = 0.125; % kPa
P.sigma_static      = 0.01;  % kPa
P.beta_diff         = 0.02;  % kPa
P.sigma_diff        = 0.002; % kPa

P.sigma_magnet      = 0.3*pi/180; % deg * pi / 180 = radians
P.beta_magnet       = 1*pi/180;   % deg * pi / 180 = radians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params for Aersonade UAV
%physical parameters of airframe
P.mass = 25.0;
P.Jx   = 0.8244;
P.Jy   = 1.135;
P.Jz   = 1.759;
P.Jxz  = .1204;
% aerodynamic coefficients
P.S_wing        = 0.55;
P.b             = 2.8956;
P.c             = 0.18994;
P.S_prop        = 0.2027;
P.rho           = 1.2682;
P.k_motor       = 80;
P.k_T_P         = 0;
P.k_Omega       = 0;
P.e             = 0.9;

P.C_L_0         = 0.28;
P.C_L_alpha     = 3.45;
P.C_L_q         = 0.0;
P.C_L_delta_e   = -0.36;
P.C_D_0         = 0.03;
P.C_D_alpha     = 0.30;
P.C_D_p         = 0.0437;
P.C_D_q         = 0.0;
P.C_D_delta_e   = 0.0;
P.C_m_0         = -0.02338;
P.C_m_alpha     = -0.38;
P.C_m_q         = -3.6;
P.C_m_delta_e   = -0.5;
P.C_Y_0         = 0.0;
P.C_Y_beta      = -0.98;
P.C_Y_p         = 0.0;
P.C_Y_r         = 0.0;
P.C_Y_delta_a   = 0.0;
P.C_Y_delta_r   = -0.17;
P.C_ell_0       = 0.0;
P.C_ell_beta    = -0.12;
P.C_ell_p       = -0.26;
P.C_ell_r       = 0.14;
P.C_ell_delta_a = 0.08;
P.C_ell_delta_r = 0.105;
P.C_n_0         = 0.0;
P.C_n_beta      = 0.25;
P.C_n_p         = 0.022;
P.C_n_r         = -0.35;
P.C_n_delta_a   = 0.06;
P.C_n_delta_r   = -0.032;
P.C_prop        = 1.0;
P.M             = 50;
P.epsilon       = 0.1592;
P.alpha0        = 0.4712;

% wind parameters
P.wind_n = 0;%3;
P.wind_e = 0;%2;
P.wind_d = 0;
P.L_u = 200;
P.L_v = 200;
P.L_w = 50;
P.sigma_u = 1.06; 
P.sigma_v = 1.06;
P.sigma_w = .7;


% compute trim conditions using 'mavsim_chap5_trim.slx'
% initial airspeed
P.Va0    = 35;
P.gamma0 = 0*pi/180;        % desired flight path angle (radians)
P.R0     = inf;      % desired radius (m) - use (+) for right handed orbit, 

% autopilot sample rate
P.Ts = 0.01;
P.Ts_gps = 1.0;

% first cut at initial conditions
P.pn0    = 0;  % initial North position
P.pe0    = 0;  % initial East position
P.pd0    = 0;  % initial Down position (negative altitude)
P.u0     = P.Va0; % initial velocity along body x-axis
P.v0     = 0;  % initial velocity along body y-axis
P.w0     = 0;  % initial velocity along body z-axis
P.phi0   = 0;  % initial roll angle
P.theta0 = 0;  % initial pitch angle
P.psi0   = 0;  % initial yaw angle
P.p0     = 0;  % initial body frame roll rate
P.q0     = 0;  % initial body frame pitch rate
P.r0     = 0;  % initial body frame yaw rate

                    %                          (-) for left handed orbit

%% run trim commands
[x_trim, u_trim]=compute_trim('mavsim_trim',P.Va0,P.gamma0,P.R0);
P.u_trim = u_trim;
P.x_trim = x_trim;

% set initial conditions to trim conditions
% initial conditions
P.pn0    = 0;  % initial North position
P.pe0    = 0;  % initial East position
P.pd0    = 0;  % initial Down position (negative altitude)
P.u0     = x_trim(4);  % initial velocity along body x-axis
P.v0     = x_trim(5);  % initial velocity along body y-axis
P.w0     = x_trim(6);  % initial velocity along body z-axis
P.phi0   = x_trim(7);  % initial roll angle
P.theta0 = x_trim(8);  % initial pitch angle
P.psi0   = x_trim(9);  % initial yaw angle
P.p0     = x_trim(10);  % initial body frame roll rate
P.q0     = x_trim(11);  % initial body frame pitch rate
P.r0     = x_trim(12);  % initial body frame yaw rate

%% compute different transfer functions
[T_phi_delta_a,T_chi_phi,T_theta_delta_e,T_h_theta,T_h_Va,T_Va_delta_t,T_Va_theta,T_v_delta_r]...
    = compute_tf_model(x_trim,u_trim,P);

% linearize the equations of motion around trim conditions
[A_lon, B_lon, A_lat, B_lat] = compute_ss_model('mavsim_trim',x_trim,u_trim);

% Compute the eigen values of A_lon and A_lat
e_lat = eig(A_lat);
e_lon = eig(A_lon);

% Find the natural frequency and damping ratio from the eigenvalues for
% longitunal motion
wn_phugoid = norm(e_lon(2));
wn_short_period = norm(e_lon(4));
zeta_phugoid = (-real(e_lon(2)))/wn_phugoid; 
zeta_short_period = (-real(e_lon(4)))/wn_short_period; 

if(wn_short_period < wn_phugoid)
    temp = wn_phugoid;
    wn_phugoid = wn_short_period;
    wn_short_period = temp;
    
    temp = zeta_phugoid;
    zeta_phugoid = zeta_short_period;
    zeta_short_period = temp;
end

% Find the natural frequency and damping ratio from the eigenvalues for
% lateral motion
for i = 2:5
   if ((imag(e_lat(i)) == 0.0) && (real(e_lat(i) > 0)) )  
       wn_spiral_divergence_mode = norm(e_lat(i));
       zeta_spiral_divergence_mode = 1;
   end
   
   if ((imag(e_lat(i)) == 0.0) && (real(e_lat(i) < 0)) )  
       wn_roll_mode = norm(e_lat(i));
       zeta_roll_mode = 1;
   end
   
   if ((imag(e_lat(i)) ~= 0.0))  
       wn_dutch_roll_mode = norm(e_lat(i));
       zeta_dutch_roll_mode = (-real(e_lat(i)))/wn_dutch_roll_mode;
   end
end

% constants for computing gains
P.delta_a_max = 45*pi/180;
P.phi_max = 30*pi/180;
P.delta_e_max = 45*pi/180;
P.theta_max = 15*pi/180;
P.delta_t_max = 1.0;
% compute gains and store them in P
[P.kp_phi, P.kd_phi, P.ki_phi, P.kp_chi, P.kd_chi, P.ki_chi, P.kp_beta, P.kd_beta, P.ki_beta, P.kp_theta, P.kd_theta, P.ki_theta, P.kp_h, P.kd_h, P.ki_h, P.kp_v2, P.kd_v2, P.ki_v2, P.kp_v, P.kd_v, P.ki_v] = compute_gains(T_phi_delta_a,T_chi_phi,T_theta_delta_e,T_h_theta,T_h_Va,T_Va_delta_t,T_Va_theta,T_v_delta_r,P);

%altitude parameters
% altitude parameters and gains
P.altitude_take_off_zone = 20;
P.altitude_hold_zone = 10;
P.theta_takeoff = 15*pi/180;
P.delta_t_climb = 0.65;