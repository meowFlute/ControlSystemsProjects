% estimate_states
% - estimate the MAV states using gyros, accels, pressure sensors, and GPS.

function xhat = estimate_states(uu, P)

    % sensor inputs
    y_gyro_x      = uu(1);
    y_gyro_y      = uu(2);
    y_gyro_z      = uu(3);
    y_accel_x     = uu(4);
    y_accel_y     = uu(5);
    y_accel_z     = uu(6);
    y_static_pres = uu(7);
    y_diff_pres   = uu(8);
    y_gps_n       = uu(9);
    y_gps_e       = uu(10);
    y_gps_h       = uu(11);
    y_gps_Vg      = uu(12);
    y_gps_course  = uu(13);
    t             = uu(14);
    
    %% GYRO LOW PASS FILTER FOR ANGULAR RATES
    persistent phat_dl;
    persistent qhat_dl;
    persistent rhat_dl;
    if isempty(phat_dl),
        phat_dl = 0;
        qhat_dl = 0;
        rhat_dl = 0;
    end
    gyro_cutoff_frequency = 50; % Hz
    roll_cutoff_frequency = 50; % Hz
    gyro_alpha_LPF = exp(-gyro_cutoff_frequency*P.Ts);
    gyro_roll_LPF = exp(-roll_cutoff_frequency*P.Ts);
    
    phat = gyro_roll_LPF*phat_dl + (1 - gyro_roll_LPF)*y_gyro_x;
    qhat = gyro_alpha_LPF*qhat_dl + (1 - gyro_alpha_LPF)*y_gyro_y;
    rhat = gyro_alpha_LPF*rhat_dl + (1 - gyro_alpha_LPF)*y_gyro_z;
    
    phat_dl = phat;
    qhat_dl = qhat;
    rhat_dl = rhat;
    
    %% ALTITUDE LOW PASS FILTER
    persistent hhat_dl;
    if isempty(hhat_dl),
        hhat_dl = 0;
    end
    altitude_cutoff_frequency = 50; % Hz
    altitude_alpha_LPF = exp(-altitude_cutoff_frequency*P.Ts);
    
    hhat = altitude_alpha_LPF*hhat_dl + ((1 - altitude_alpha_LPF)*((101.325-y_static_pres)*1000/(P.rho*P.gravity)));
    
    hhat_dl = hhat;
    
    %% AIRSPEED LOW PASS FILTER
    persistent Vahat_dl;
    if isempty(Vahat_dl),
        Vahat_dl = 0;
    end
    airspeed_cutoff_frequency = 50; % Hz
    airspeed_alpha_LPF = exp(-airspeed_cutoff_frequency*P.Ts);
    
    Vahat = airspeed_alpha_LPF*Vahat_dl + (1 - airspeed_alpha_LPF)*sqrt(2*y_diff_pres/P.rho);
    
    Vahat_dl = Vahat;
    
    %% ROLL AND PITCH ANGLE LOW PASS FILTER
%     phihat   = 0;
%     thetahat = 0;
%     
%     persistent filtered_accel_x_dl;
%     persistent filtered_accel_y_dl;
%     persistent filtered_accel_z_dl;
%     if isempty(filtered_accel_x_dl),
%         filtered_accel_x_dl = 0;
%         filtered_accel_y_dl = 0;
%         filtered_accel_z_dl = 0;
%     end
%     attitude_cutoff_frequency = 25; % Hz
%     attitude_alpha_LPF = exp(-attitude_cutoff_frequency*P.Ts);
%     
%     filtered_accel_x = attitude_alpha_LPF*filtered_accel_x_dl + (1 - attitude_alpha_LPF)*y_accel_x;
%     filtered_accel_y = attitude_alpha_LPF*filtered_accel_y_dl + (1 - attitude_alpha_LPF)*y_accel_y;
%     filtered_accel_z = attitude_alpha_LPF*filtered_accel_z_dl + (1 - attitude_alpha_LPF)*y_accel_z;
%     
%     filtered_accel_x_dl = filtered_accel_x;
%     filtered_accel_y_dl = filtered_accel_y;
%     filtered_accel_z_dl = filtered_accel_z;
%     
%     phihat   = atan(filtered_accel_y/filtered_accel_z);
%     thetahat = asin(filtered_accel_x/P.gravity);

    %% ROLL AND PITCH ANGLE EXTENDED KALMAN FILTER
    
    % initialize xhat as x_0 and PP
    persistent x_hat_attitude;
    persistent PP;
    if t < P.Ts,
        x_hat_attitude = [P.phi0; P.theta0;];
        PP = diag([0, 0]);
    else
        u_attitude = [phat; qhat; rhat; Vahat];
        y_accel = [y_accel_x; y_accel_y; y_accel_z];

        phi   = x_hat_attitude(1);
        theta = x_hat_attitude(2);
        p = u_attitude(1);
        q = u_attitude(2);
        r = u_attitude(3);
        Va = u_attitude(4);

        % the state space equations
        
        % other useful parameters    
        R = [P.sigma_gyro^2; P.sigma_gyro^2; P.sigma_gyro^2];    
        Q = diag([0,0]); %Process Noise - tune this

        % prediction step 10 times as frequent
        N = 10;
        for i = 1:N
           phi   = x_hat_attitude(1);
           theta = x_hat_attitude(2);
           f = [p + q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);...
                q*cos(phi)-r*sin(phi)];
           dfdx = [q*cos(phi)*tan(theta)-r*sin(phi)*tan(theta),   (q*sin(phi)+r*cos(phi))/cos(theta)^2;...
                  -q*sin(phi)-r*cos(phi),                          0];
              
           x_hat_attitude = x_hat_attitude + (P.Ts/N)*f;
           A = dfdx;
           PP = PP + (P.Ts/N)*(A*PP+ PP*A' + Q);
        end
        
        phi   = x_hat_attitude(1);
        theta = x_hat_attitude(2);
        dhdx = [0,                                q*Va*cos(theta)+P.gravity*cos(theta);...
               -P.gravity*cos(phi)*cos(theta),   -r*Va*sin(theta)-p*Va*cos(theta)+P.gravity*sin(phi)*sin(theta);...
                P.gravity*sin(phi)*cos(theta),    (q*Va+P.gravity*cos(phi))*sin(theta)];
            
        h = [   q*Va*sin(theta)+P.gravity*sin(theta);...
                r*Va*cos(theta)-p*Va*sin(theta)-P.gravity*cos(theta)*sin(phi);...
               -q*Va*cos(theta)-P.gravity*cos(theta)*cos(phi)];
        
        % sensor measurement update step
        for i = 1:3
           Ci = dhdx(i,:);
           Ri = R(i,:);
           Li = PP*Ci'*inv(Ri+Ci*PP*Ci');
           PP = (eye(2)-Li*Ci)*PP;
           x_hat_attitude = x_hat_attitude + Li*(y_accel(i)-h(i));
        end
    end
    
    phihat = x_hat_attitude(1);
    thetahat = x_hat_attitude(2);
    
    %% GPS LOW PASS FILTER
%     persistent pnhat_dl;
%     persistent pehat_dl;
%     persistent chihat_dl;
%     persistent Vghat_dl;
%     if isempty(pnhat_dl),
%         pnhat_dl = 0;
%         pehat_dl = 0;
%         chihat_dl = 0;
%         Vghat_dl = 0;
%     end
%     position_cutoff_frequency = 50; % Hz
%     position_alpha_LPF = exp(-position_cutoff_frequency*P.Ts);
%     
%     pnhat  = position_alpha_LPF*pnhat_dl + (1 - position_alpha_LPF)*y_gps_n;
%     pehat  = position_alpha_LPF*pehat_dl + (1 - position_alpha_LPF)*y_gps_e;
%     chihat = position_alpha_LPF*chihat_dl + (1 - position_alpha_LPF)*y_gps_course;
%     Vghat  = position_alpha_LPF*Vghat_dl + (1 - position_alpha_LPF)*y_gps_Vg;
%     
%     pnhat_dl = pnhat;
%     pehat_dl = pehat;
    
    %% GPS Extended Kalman Filter
    persistent x_hat_gps;
    persistent gps_dl;
    persistent PP_gps;
    if t < P.Ts,
        x_hat_gps = [P.pn0; P.pe0; P.Va0*cos(P.gamma0); P.psi0; P.wind_n; P.wind_e; P.psi0];
        gps_measurement = [y_gps_n; y_gps_e; y_gps_Vg; y_gps_course;];
        gps_dl = gps_measurement;
        PP_gps = diag([.1^2, .1^2, .05^2, .1^2, .4^2, .4^2, .05^2]);
%         PP_gps = diag([0, 0.1, 0.001, 0.0001, 0.01, 0.03, 0]);
    else
        pn  = x_hat_gps(1);
        pe  = x_hat_gps(2);
        Vg 	= x_hat_gps(3);
        chi = x_hat_gps(4);
        wn 	= x_hat_gps(5);
        we 	= x_hat_gps(6);
        psi = x_hat_gps(7);
        
        u_gps = [Va; qhat; rhat; phihat; thetahat];
        y_gps = [y_gps_n; y_gps_e; y_gps_Vg; y_gps_course; Va*cos(psi)+wn-Vg*cos(chi); Va*sin(psi)+we-Vg*sin(chi)];
        gps_measurement = [y_gps_n; y_gps_e; y_gps_Vg; y_gps_course;];
        
        Va    = u_gps(1);
        q     = u_gps(2);
        r     = u_gps(3);
        phi   = u_gps(4);
        theta = u_gps(5);
        
        % other useful parameters    
        R = [0.0025^2, 0.0025^2, 0.0025^2, 0.0025^2, .3^2, .3^2, .1^2]';    
        Q = diag([.1,.1,.01,0.000001,.001,.01,0]); %Process Noise - tune this

        % prediction step 10 times as frequent
        N = 10;
        for i = 1:N
            pn  = x_hat_gps(1);
            pe  = x_hat_gps(2);
            Vg 	= x_hat_gps(3);
            chi = x_hat_gps(4);
            wn 	= x_hat_gps(5);
            we 	= x_hat_gps(6);
            psi = x_hat_gps(7);
            
            psidot = q*sin(phi)/cos(theta)+r*cos(phi)/cos(theta);
            Vgdot = ((Va*cos(psi)+wn)*(-Va*psidot*sin(psi))+(Va*sin(psi)+we)*(Va*psidot*cos(psi)))/Vg;

            f = [Vg*cos(chi);...
                Vg*sin(chi);...]
                ((Va*cos(psi)+wn)*(-Va*psidot*sin(psi))+(Va*sin(psi)+we)*(Va*psidot*cos(psi)))/Vg;...
                P.g/Vg*tan(phi)*cos(chi-psi);...
                0;
                0;
                q*sin(phi)/cos(theta)+r*cos(phi)/cos(theta)];

            dfdx = [...
                0,  0,  cos(chi),                           -Vg*sin(chi),                   0,                      0,                  0;...
                0,  0,  sin(chi),                           Vg*cos(chi),                    0,                      0,                  0;...
                0,  0,  -Vgdot/Vg,                          0,                             -psidot*Va*sin(psi)/Vg,  psidot*Va*cos(psi)/Vg, (-psidot*Va*(wn*cos(psi)+we*sin(psi)))/Vg;...
                0,  0,  -P.g/Vg^2*tan(phi)*cos(chi-psi),    -P.g/Vg*tan(phi)*sin(chi-psi),  0,                      0,                  P.g/Vg*tan(phi)*sin(chi-psi);...
                0,  0,  0,                                  0,                              0,                      0,                  0;...
                0,  0,  0,                                  0,                              0,                      0,                  0;...
                0,  0,  0,                                  0,                              0,                      0,                  0];
                
           x_hat_gps = x_hat_gps + (P.Ts/N)*f;
           A = dfdx;
           PP_gps = PP_gps + (P.Ts/N)*(A*PP_gps+ PP_gps*A' + Q);
        end
        pn  = x_hat_gps(1);
        pe  = x_hat_gps(2);
        Vg 	= x_hat_gps(3);
        chi = x_hat_gps(4);
        wn 	= x_hat_gps(5);
        we 	= x_hat_gps(6);
        psi = x_hat_gps(7);
        
        dhdx = [...
                1,  0,  0,          0,              0,  0,  0;...
                0,  1,  0,          0,              0,  0,  0;...
                0,  0,  1,          0,              0,  0,  0;...
                0,  0,  0,          1,              0,  0,  0;...
                0,  0,  -cos(chi),  Vg*sin(chi),    1,  0,  -Va*sin(psi);...
                0,  0,  -sin(chi),  -Vg*cos(chi),   0,  1,  Va*cos(psi)];
            
        h = [   pn; pe; Vg; chi; Va*cos(psi)+wn-Vg*cos(chi); Va*sin(psi) + we - Vg*sin(chi) ];
        % sensor measurement update step when you recieve data
        if(any(gps_dl ~= gps_measurement)),
            for i = 1:6
               Ci = dhdx(i,:);
               Ri = R(i,:);
               Li = PP_gps*Ci'*inv(Ri+Ci*PP_gps*Ci');
               PP_gps = (eye(7)-Li*Ci)*PP_gps;
               x_hat_gps = x_hat_gps + Li*(y_gps(i)-h(i));
            end
            gps_dl = gps_measurement;
        end
    end
    
    pnhat  = x_hat_gps(1);
    pehat  = x_hat_gps(2);
    Vghat 	= x_hat_gps(3);
    chihat = x_hat_gps(4);
    wnhat 	= x_hat_gps(5);
    wehat 	= x_hat_gps(6);
    psihat = x_hat_gps(7);
  
    %% OUTPUT
    % not estimating these states 
    alphahat = 0;
    betahat  = 0;
    bxhat    = 0;
    byhat    = 0;
    bzhat    = 0;
    
      xhat = [...
        pnhat;...
        pehat;...
        hhat;...
        Vahat;...
        alphahat;...
        betahat;...
        phihat;...
        thetahat;...
        chihat;...
        phat;...
        qhat;...
        rhat;...
        Vghat;...
        wnhat;...
        wehat;...
        psihat;...
        bxhat;...
        byhat;...
        bzhat;...
        ];
end
