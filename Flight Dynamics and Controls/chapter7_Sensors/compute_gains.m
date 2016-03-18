function [kp_phi, kd_phi, ki_phi, kp_chi, kd_chi, ki_chi, kp_beta, kd_beta, ki_beta, kp_theta, kd_theta, ki_theta, kp_h, kd_h, ki_h, kp_v2, kd_v2, ki_v2, kp_v, kd_v, ki_v] = compute_gains(T_phi_delta_a,T_chi_phi,T_theta_delta_e,T_h_theta,T_h_Va,T_Va_delta_t,T_Va_theta,T_v_delta_r,P)
    %% Problem 6.1 - Roll Hold Loop Gains
    %  Create a matlab script that computes the gains for the roll attitude
    %  hold loop. 
    
                    delta_a_max = P.delta_a_max;
                    phi_max     = P.phi_max;
                    
    % Using the transfer function of the open loop we know a_phi 1 and 2
    [numerator, denomenator] = tfdata(T_phi_delta_a, 'v');
    a_phi_2 = numerator(3);
    a_phi_1 = denomenator(2);
    
    % Section 6.3 then details how to use that info to tune the roll loop
    kp_phi = (delta_a_max/phi_max)*sign(a_phi_2);
    
    zeta_phi = 0.707;
    wn_phi   = sqrt(abs(a_phi_2*(delta_a_max/phi_max)));
    
    kd_phi = (2*zeta_phi*wn_phi - a_phi_1)/a_phi_2;
    ki_phi = 0;
    
    %% Problem 6.2 - Course Hold Loop Gains
    %  tune these parameters to get 
    
    % here are our tuning parameters
    W_chi = 10;             % Pick something bigger than 5
    wn_chi = wn_phi/W_chi;  % doing so picks our natural frequency
    zeta_chi = 0.707;       % this is just a tuning parameter for damping
    Vg = P.Va0*cos(P.gamma0);
    
    % here are the gains
    kp_chi = 2*zeta_chi*wn_chi*Vg/P.gravity;
    kd_chi = 0;
    ki_chi = (wn_chi^2)*Vg/P.gravity;
    
    %% Problem 6.3 - Sideslip hold
    [numerator, denomenator] = tfdata(T_v_delta_r, 'v');
    a_beta_2 = numerator(2);
    a_beta_1 = denomenator(2);
    
    % tuning parameters
    delta_r_max = 45*pi/180;
    beta_max = 5*pi/180;
    zeta_beta = 0.707;
    
    % gains
    kp_beta = (delta_r_max/beta_max)*sign(a_beta_2);
    kd_beta = 0;
    ki_beta = (1/a_beta_2)*(((a_beta_1 + a_beta_2*kp_beta)/(2*zeta_beta))^2);
    
    %% Problem 6.4 - Pitch attitude hold loop
    
    % tuning parameters
    delta_e_max = P.delta_e_max;
    theta_max = P.theta_max;
    zeta_theta = 0.707;
    
    [numerator, denomenator] = tfdata(T_theta_delta_e, 'v');
    a_theta_3 = numerator(3);
    a_theta_1 = denomenator(2);
    a_theta_2 = denomenator(3);
    
    % gains
    kp_theta = (delta_e_max/theta_max)*sign(a_theta_3);
    
    wn_theta = sqrt(a_phi_2 + ((delta_e_max/theta_max)*abs(a_theta_3)));
    
    kd_theta = (2*zeta_theta*wn_theta - a_theta_1)/a_theta_3;
    ki_theta = 0;
    
    % we'll need to know the DC gain later since it isn't necessarily 1
    K_theta_DC = kp_theta*a_theta_3/(a_theta_2 + kp_theta*a_theta_3);
    
    %% Problem 6.5
    
    % tuning parameters
    Wh = 30; % usually between 5 and 15
    wn_h = wn_theta/Wh;
    zeta_h = 1.5;
    
    %gains
    kp_h = 2*zeta_h*wn_h/(K_theta_DC*P.Va0);
    kd_h = 0;
    ki_h = (wn_h^2)/(K_theta_DC*P.Va0);
    
    %% Problem 6.6 airspeed from pitch
    
    % tuning parameters
    WV_2 = 10;
    wn_v_2 = wn_theta/WV_2;
    zeta_v_2 = 0.707;
    
    % constants from the transfer function model
    [numerator, denomenator] = tfdata(T_Va_theta, 'v');
    a_v_1 = denomenator(2);
    a_v_3 = -numerator(2);
    
    % gains
    kp_v2 = (a_v_1 - 2*zeta_v_2*wn_v_2)/(K_theta_DC*P.gravity);
    kd_v2 =  0;
    ki_v2 =  -(wn_v_2^2)/(K_theta_DC*P.gravity);
    
    %% Problem 6.7 airspeed from throttle
    [numerator, denomenator] = tfdata(T_Va_delta_t, 'v');
    a_v_2 = numerator(2);
    
    %  tuning parameters
    WV = 10;
    wn_v = wn_theta/WV;
    zeta_v = 0.707;
    
    %  gains
    kp_v = (2*zeta_v*wn_v - a_v_1)/(a_v_2);
    kd_v = 0;
    ki_v = (wn_v^2)/(a_v_2);
    
end