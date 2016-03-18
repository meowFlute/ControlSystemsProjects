% forces_moments.m
%   Computes the forces and moments acting on the airframe. 
%
%   Output is
%       F     - forces
%       M     - moments
%       Va    - airspeed
%       alpha - angle of attack
%       beta  - sideslip angle
%       wind  - wind vector in the inertial frame
%

function out = forces_moments(x, delta, wind, P)

    % relabel the inputs
    pn      = x(1);
    pe      = x(2);
    pd      = x(3);
    u       = x(4);
    v       = x(5);
    w       = x(6);
    phi     = x(7);
    theta   = x(8);
    psi     = x(9);
    p       = x(10);
    q       = x(11);
    r       = x(12);
    delta_e = delta(1);
    delta_a = delta(2);
    delta_r = delta(3);
    delta_t = delta(4);
    w_ns    = wind(1); % steady wind - North
    w_es    = wind(2); % steady wind - East
    w_ds    = wind(3); % steady wind - Down
    u_wg    = wind(4); % gust along body x-axis
    v_wg    = wind(5); % gust along body y-axis    
    w_wg    = wind(6); % gust along body z-axis
    
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
    Rv_b = rot3*rot2*rot1;  % from velocity to body frame
    
    % compute total wind in body frame
    Vb_w = Rv_b*[w_ns; w_es; w_ds;] + [u_wg; v_wg; w_wg;];
    
    % compute wind data in NED ...why?
    Vv_w = Rv_b'*Vb_w;
    w_n = Vv_w(1);
    w_e = Vv_w(2);
    w_d = Vv_w(3);
    
    % compute body-frame components of airspeed
    Vb_a = [u; v; w;] - Vb_w;
    u_r = Vb_a(1);
    v_r = Vb_a(2);
    w_r = Vb_a(3);
    
    % compute air data
    Va = sqrt(u_r^2 + v_r^2 + w_r^2);
    alpha = atan(w_r/u_r);
    beta = asin(v_r/Va);
    
    %define all of the Cx and Cz functions
    C_L_alpha = P.C_L_0 + P.C_L_alpha*alpha;
    C_D_alpha = P.C_D_0 + P.C_D_alpha*alpha;
    C_X_alpha = -C_D_alpha*cos(alpha) + C_L_alpha*sin(alpha);
    C_X_q_alpha = -P.C_D_q*cos(alpha) + P.C_L_q*sin(alpha);
    C_X_delta_e_alpha = -P.C_D_delta_e*cos(alpha) + P.C_L_delta_e*sin(alpha);
    C_Z_alpha = -C_D_alpha*sin(alpha) - C_L_alpha*cos(alpha);
    C_Z_q_alpha = -P.C_D_q*sin(alpha) - P.C_L_q*cos(alpha);
    C_Z_delta_e_alpha = -P.C_D_delta_e*sin(alpha) - P.C_L_delta_e*cos(alpha);
    
    % compute external forces and torques on aircraft
            % double check on this matrix complete
    Force = [  -P.mass*P.gravity*sin(theta);... 
                P.mass*P.gravity*cos(theta)*sin(phi);...
                P.mass*P.gravity*cos(theta)*cos(phi);]...
         ...% double check on this matrix complete
          + ((1/2)*P.rho*Va^2*P.S_wing*...
            [   C_X_alpha + C_X_q_alpha*(P.c/(2*Va))*q + C_X_delta_e_alpha*delta_e;...
                P.C_Y_0 + P.C_Y_beta*beta + P.C_Y_p*(P.b/(2*Va))*p + P.C_Y_r*(P.b/(2*Va))*r + P.C_Y_delta_a*delta_a + P.C_Y_delta_r*delta_r;...
                C_Z_alpha + C_Z_q_alpha*(P.c/(2*Va))*q + C_Z_delta_e_alpha*delta_e])...
         ...% double check on this matrix complete
          + ((1/2)*P.rho*P.S_prop*P.C_prop*...
            [   (P.k_motor*delta_t)^2 - Va^2;...
                 0;...
                 0;]);
            % double check on this matrix 
    Torque = (1/2)*P.rho*Va^2*P.S_wing*...
            [   P.b*(P.C_ell_0 + P.C_ell_beta*beta + P.C_ell_p*(P.b/(2*Va))*p + P.C_ell_r*(P.b/(2*Va))*r + P.C_ell_delta_a*delta_a + P.C_ell_delta_r*delta_r);...
                P.c*(P.C_m_0 + P.C_m_alpha*alpha + P.C_m_q*(P.c/(2*Va))*q + P.C_m_delta_e*delta_e);...
                P.b*(P.C_n_0 + P.C_n_beta*beta + P.C_n_p*(P.b/(2*Va))*p + P.C_n_r*(P.b/(2*Va))*r + P.C_n_delta_a*delta_a + P.C_n_delta_r*delta_r);...
                ]...
            + [-P.k_T_P*(P.k_Omega*delta_t)^2; 0; 0;];
                
    out = [Force; Torque; Va; alpha; beta; w_n; w_e; w_d];
end



