loadDronePoints;
P.gravity = 9.8;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params for Aersonade UAV
%physical parameters of airframe
P.mass = 25;
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
P.wind_n = 3;
P.wind_e = 2;
P.wind_d = 0;
P.L_u = 200;
P.L_v = 200;
P.L_w = 50;
P.sigma_u = 1.06; 
P.sigma_v = 1.06;
P.sigma_w = .7;

% compute trim conditions using 'mavsim_chap5_trim.slx'
% initial airspeed
P.Va0 = 35;        % m/s (~85 mph)
gamma = 15*pi/180;  % desired flight path angle (radians)
R     = Inf;       % desired radius (m) - use (+) for right handed orbit, 
                   %                          (-) for left handed orbit
                   
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

% run trim commands
[x_trim, u_trim, y_trim]=compute_trim('mavsim_trim',P.Va0,gamma,R);
P.u_trim = u_trim;
P.x_trim = x_trim;
P.y_trim = y_trim;

% autopilot sample rate
P.Ts = 0.01;

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

% compute different transfer functions
[T_phi_delta_a,T_chi_phi,T_theta_delta_e,T_h_theta,T_h_Va,T_Va_delta_t,T_Va_theta,T_v_delta_r]...
    = compute_tf_model(x_trim,u_trim,y_trim,P);

% linearize the equations of motion around trim conditions
[A_lon, B_lon, A_lat, B_lat] = compute_ss_model('mavsim_trim',x_trim,u_trim);

% PID constants

% roll loop
[num, den] = tfdata(T_phi_delta_a,'v');
a_phi_2 = num(3);
a_phi_1 = den(2); 
P.delta_a_max = 45 * pi()/180;
P.phi_max = 30 * pi()/180;
P.kp_roll = P.delta_a_max / P.phi_max * sign(a_phi_2); % kp_roll
roll_zeta = 1.5;
roll_wn = sqrt(abs(a_phi_2) * P.delta_a_max / P.phi_max);
P.kd_roll = (2*roll_zeta*roll_wn-a_phi_1) / a_phi_2; % kd_roll

% course loop
[num, den] = tfdata(T_chi_phi,'v');
gVg = num(end);
course_wn = 1/(15) * roll_wn;
P.ki_course = course_wn^2 / gVg;
course_zeta = 0.707;
P.kp_course = 2*course_zeta*course_wn / gVg;

% pitch loop
%
[num, den] = tfdata(T_theta_delta_e,'v');
a_theta_1 = den(2);
a_theta_2 = den(3);
a_theta_3 = num(end);

P.delta_e_max = 45 * pi()/180;
P.theta_max = 30 * pi()/180;
P.kp_pitch = P.delta_e_max / P.theta_max * sign(a_theta_3);
pitch_wn = sqrt(a_theta_2 + P.kp_pitch * a_theta_3);
pitch_zeta = 0.707;
P.kd_pitch = (2 * pitch_zeta * pitch_wn - a_theta_1) / a_theta_3;
K_theta_DC = (P.kp_pitch * a_theta_3) / (a_theta_2 + P.kp_pitch * a_theta_3);
P.K_theta_DC = K_theta_DC;

% airspeed w throttle hold loop
P.delta_t_max = 1;
Va_trim = P.y_trim(1);
alpha_trim = P.y_trim(2);
delta_e_trim = P.u_trim(1);
delta_t_trim = P.u_trim(4);
[num, den] = tfdata(T_Va_delta_t,'v');
a_v_1 = den(end);
a_v_2 = num(end);
[num, den] = tfdata(T_Va_theta,'v');
a_v_3 = num(end);

wn_v = 5; % tuning parameter
zeta_v = 0.707; % tuning parameter
P.ki_v = wn_v^2 / a_v_2;
P.kp_v = (2 * zeta_v * wn_v - a_v_1) / a_v_2;

% airspeed w pitch hold loop
wn_v2 = 1/(10) * pitch_wn;
P.ki_v2 = -wn_v2^2 / K_theta_DC / P.gravity;
zeta_v2 = 0.707;
P.kp_v2 = (a_v_1 - 2*zeta_v2*wn_v2) / K_theta_DC / P.gravity;

% altitude w pitch hold loop
wn_h = 1 / (15) * pitch_wn;
P.ki_h = wn_h^2 / (K_theta_DC * P.Va0);
zeta_h = 0.707;
P.kp_h = (2*zeta_h*wn_h) / (K_theta_DC * P.Va0);

% take off/cruise parameters
P.altitude_take_off_zone = 30;
P.altitude_hold_zone = 10;

%----------------SENSORS----------------
P.sigma_gyro = 0.13 * pi()/180;
P.sigma_accel = 0.0025;
P.beta_abs_pres = 0.125;
P.sigma_abs_pres = 0.01;
P.beta_diff_pres = 0.02;
P.sigma_diff_pres = 0.002;
P.k_gps = 1/1100;
P.Ts_gps = 1.0;
P.sigma_n_gps = 0.21;
P.sigma_e_gps = 0.21;
P.sigma_h_gps = 0.4;

P.bias_gyro_x = 0;
P.bias_gyro_y = 0;
P.bias_gyro_z = 0;


% need to add the following to your parameter file:

% chapter 11 - path manager
% number of waypoints in data structure
h0 = P.pd0;
P.size_waypoint_array = 10;
P.R_min = P.Va0^2/P.gravity/tan(P.phi_max);

% create random city map
city_width      = 2000;  % the city is of size (width)x(width)
building_height = 300;   % maximum height of buildings
%building_height = 1;   % maximum height of buildings (for camera)
num_blocks      = 5;    % number of blocks in city
street_width    = .8;   % percent of block that is street.
P.pd0           = -h0;  % initial height of MAV
P.map = createWorld(city_width, building_height, num_blocks, street_width);

