% select gains for roll loop
    % get transfer function data for delta_a to phi
    [num,den]=tfdata(T_phi_delta_a,'v');
    a_phi2 = num(3);
    a_phi1 = den(2);
    % maximum possible aileron command
    P.delta_a_max = 30*pi/180;
    % Roll command when delta_a_max is achieved
    P.phi_max = 45*pi/180;
    % pick natural frequency to achieve delta_a_max for step of phi_max
    zeta_roll = 0.707;
    wn_roll = sqrt(abs(a_phi2)*P.delta_a_max/P.phi_max);%sqrt(a_phi2*delta_a_max*sqrt(1-zeta_roll^2)/phi_max);
    
    % set control gains based on zeta and wn
    P.roll_kp = P.delta_a_max/P.phi_max*sign(a_phi2);
    P.roll_kd = (2*zeta_roll*wn_roll-a_phi1)/a_phi2;
    P.roll_ki = 0; %We don't want an integrator on this inner loop
    
   
    
% select gains for course loop
    %design parameters
    
    [num,den] = tfdata(T_chi_phi,'v');
    
    gVg = num(end);
    
    W_course = 10; %Bandwitdh separation
    zeta_course = .707;  %damping ratio
    
   wn_course = 1/W_course*wn_roll;
   P.course_kp = 2*zeta_course*wn_course/gVg;
   P.course_ki = wn_course^2/gVg;
   P.course_kd = 0;  % No differentiator for this loop
   
% select gains for sideslip hold
    P.beta_kp = 0;
    P.beta_ki = 0;
    P.beta_kd = 0;

   
% select gains for the pitch loop

    [num,den] = tfdata(T_theta_delta_e, 'v');
    a_theta3 = num(3);
    a_theta2 = den(3);
    a_theta1 = den(2);
    
    P.delta_e_max = 45 * pi/180;
    P.theta_max = 15 * pi/180;
    
    wn_pitch = sqrt(a_theta2+(P.delta_e_max/P.theta_max)*abs(a_theta3));
    zeta_pitch = .707;
    
   P.pitch_kp = P.delta_e_max/P.theta_max*sign(a_theta3);
   P.pitch_kd = (2*zeta_pitch*wn_pitch-a_theta1)/a_theta3;
   P.pitch_ki = 0.0;  %No integrator for this loop
   P.K_theta_DC = P.pitch_kp*a_theta3/(a_theta2+P.pitch_kp*a_theta3);

% select gains for altitude loop

    W_altitude = 10;
    zeta_altitude = .707;
    
    wn_altitude = 1/W_altitude*wn_pitch;
    
   P.altitude_kp = (2*zeta_altitude*wn_altitude)/(P.K_theta_DC*P.Va);
   P.altitude_ki = wn_altitude^2/(P.K_theta_DC*P.Va);
   P.altitude_kd = 0;
 
% airspeed hold using pitch

    [num,den] = tfdata(T_Va_theta,'v');
    a_v1 = den(end);
    
    Wv2 = 10;
    wn_v2 = wn_pitch/Wv2;
    zeta_v2 = .707;
    
    
   P.airspeed_pitch_kp = (a_v1-2*zeta_v2*wn_v2)/(P.K_theta_DC*P.g);
   P.airspeed_pitch_ki = -wn_v2^2/(P.K_theta_DC*P.g);
   P.airspeed_pitch_kd = 0;
 
% airspeed hold using throttle
[num,den] = tfdata(T_Va_delta_t,'v');
a_v2 = num(end);

    %design parameters
        wn_v = .5;
        zeta_v = .707;


   P.airspeed_throttle_kp = (2*zeta_v*wn_v-a_v1)/a_v2;
   P.airspeed_throttle_ki = wn_v^2/a_v2;
   P.airspeed_throttle_kd = 0;
 


