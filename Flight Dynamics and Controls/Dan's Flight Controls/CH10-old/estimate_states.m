% estimate_states
%   - estimate the MAV states using gyros, accels, pressure sensors, and
%   GPS.
%
% Outputs are:
%   pnhat    - estimated North position, 
%   pehat    - estimated East position, 
%   hhat     - estimated altitude, 
%   Vahat    - estimated airspeed, 
%   alphahat - estimated angle of attack
%   betahat  - estimated sideslip angle
%   phihat   - estimated roll angle, 
%   thetahat - estimated pitch angel, 
%   chihat   - estimated course, 
%   phat     - estimated roll rate, 
%   qhat     - estimated pitch rate, 
%   rhat     - estimated yaw rate,
%   Vghat    - estimated ground speed, 
%   wnhat    - estimate of North wind, 
%   wehat    - estimate of East wind
%   psihat   - estimate of heading angle
% 
% 
% Modified:  3/15/2010 - RB
%            5/18/2010 - RB
%

function xhat = estimate_states(uu, P)

   % rename inputs
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
   
    a = 150;
    agps = 200;
    
   LPF = @(a, ts, x, x1) exp(-a*ts)*x1 + (1-exp(-a*ts))*x;
   persistent gyrox1 gyroy1 gyroz1
   persistent staticpres1 diffpres1
   persistent accelx1 accely1 accelz1
   persistent wnhat wehat;
   persistent gpsv1;
   persistent pnhat pehat Vghat chihat
   persistent phihat thetahat
   if(t == 0)
      gyrox1 = y_gyro_x;
      gyroy1 = y_gyro_y;
      gyroz1 = y_gyro_z;
      staticpres1 = y_static_pres;
      diffpres1 = y_diff_pres;
      accelx1 = y_accel_x;
      accely1 = y_accel_y;
      accelz1 = y_accel_z;
      gpsv1 = y_gps_Vg;
      wnhat = 0;
      wehat = 0;
      pnhat = 0;
      pehat = 0;
      chihat = 0;
      phihat = 0;
      thetahat = 0;
   end
   
    phat = LPF(a, P.Ts, y_gyro_x, gyrox1);
    qhat = LPF(a, P.Ts, y_gyro_y, gyroy1);
    rhat = LPF(a, P.Ts, y_gyro_z, gyroz1);
    hhat = LPF(a, P.Ts, y_static_pres, staticpres1)/P.rho/P.g;
    Vahat = sqrt(2/P.rho*LPF(a, P.Ts, y_diff_pres, diffpres1));

%      phihat = atan(LPF(a, P.Ts, y_accel_y, accely1)/LPF(a, P.Ts, y_accel_z, accelz1));
%      thetahat = asin(LPF(a, P.Ts, y_accel_x, accelx1)/P.g);

    Vghat = LPF(agps, P.Ts_gps, gpsv1, y_gps_Vg);
    psihat = chihat;
    
    
    %Kalman Filter: Attitude
    xhat_att = [phihat; thetahat];
    u_att = [phat; qhat; rhat; Vahat];
    y_accel = [y_accel_x; y_accel_y; y_accel_z];
    xhat_att = EKF1(xhat_att, u_att, y_accel, P);
    
    phihat = xhat_att(1);
    thetahat = xhat_att(2);
    
    %Kalman Filter: GPS
    xhat_gps = [pnhat; pehat; Vghat; chihat; wnhat; wehat; psihat];
    y_gps = [y_gps_n; y_gps_e; y_gps_Vg; y_gps_course];
    u_gps = [u_att; phihat; thetahat];
    xhat_gps = EKF2(xhat_gps, u_gps, y_gps, P, t);
    
    pnhat = xhat_gps(1);
    pehat = xhat_gps(2);
    Vghat = xhat_gps(3);
    chihat = xhat_gps(4);
    wnhat = xhat_gps(5);
    wehat = xhat_gps(6);
    psihat = xhat_gps(7);
    
%     pnhat = 0;
%     pehat = 0;
%     Vghat = 0;
%     chihat = 0;
%     wnhat = 0;
%     wehat = 0;
%     psihat = 0;
   
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

function xhat = EKF1(xhat, u, y_accel, P)
    persistent PP
    if(isempty(PP))
        PP = diag([.02^2, .02^2]);
    end
    
    phi = xhat(1);
    theta = xhat(2);

    p = u(1);
    q = u(2);
    r = u(3);
    Va = u(4);

    N = 10;
    
    f = [p + q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);...
        q*cos(phi)-r*sin(phi)];
    dfdx = [...
            q*cos(phi)*tan(theta)-r*sin(phi)*tan(theta),   (q*sin(phi)-r*cos(phi))/cos(theta)^2;...
            -q*sin(phi)-r*cos(phi),                         0];
        
    Q = diag([.0001,.0001]); %Process Noise - tune this

    for i = 1:N
       xhat = xhat + (P.Ts/N)*f;
       A = dfdx;
       PP = PP + (P.Ts/N)*(A*PP+PP*A'+Q);
    end
    
    dhdx = [...
            0,                          q*Va*cos(theta)+P.g*cos(theta);...
            -P.g*cos(phi)*cos(theta),   -r*Va*sin(theta)-p*Va*cos(theta)+P.g*sin(phi)*sin(theta);...
            P.g*sin(phi)*cos(theta),    (q*Va+P.g*cos(phi))*sin(theta)];
        
    h = [...
        q*Va*sin(theta)+P.g*sin(theta);...
        r*Va*cos(theta)-p*Va*sin(theta)-P.g*cos(theta)*sin(phi);...
        -q*Va*cos(theta)-P.g*cos(theta)*cos(phi)];
    
    R = [P.sigma_gyro^2; P.sigma_gyro^2; P.sigma_gyro^2];
    
    for i = 1:3
       Ci = dhdx(i,:);
       Ri = R(i,:);
       Li = PP*Ci'/(Ri+Ci*PP*Ci');
       PP = (eye(2)-Li*Ci)*PP;
       xhat = xhat + Li*(y_accel(i)-h(i));
    end
    
end

function xhat = EKF2(xhat, u, gps, P, t)

    persistent gps_d1 PP;
    if(t == 0)
       gps_d1 = [0; 0; 0; 0]; 
       PP = diag([.1^2, .1^2, .05^2, .1^2, .4^2, .4^2, .05^2]);
    end
    N = 10;
    y_gps_n = gps(1);
    y_gps_e = gps(2);
    y_gps_Vg = gps(3);
    y_gps_course = gps(4);
    
    pn = xhat(1);
    pe = xhat(2);
    Vg = xhat(3);
    chi = xhat(4);
    wn = xhat(5);
    we = xhat(6);
    psi = xhat(7);
    
    p = u(1);
    q = u(2);
    r = u(3);
    Va = u(4);
    phi = u(5);
    theta = u(6);

    
    for i= 1:N
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
            0,  0,  -Vgdot/Vg,                          0,                              -psidot*Va*sin(psi),    psidot*Va*cos(psi), (-psidot*Va*(wn*cos(psi)+we*sin(psi)))/Vg;...
            0,  0,  -P.g/Vg^2*tan(phi)*cos(chi-psi),    -P.g/Vg*tan(phi)*sin(chi-psi),  0,                      0,                  P.g/Vg*tan(phi)*sin(chi-psi);...
            0,  0,  0,                                  0,                              0,                      0,                  0;...
            0,  0,  0,                                  0,                              0,                      0,                  0;...
            0,  0,  0,                                  0,                              0,                      0,                  0];
        
        Q = diag([.2,.2,.2,.2,.2,.2,.2]);
        xhat = xhat+(P.Ts/N)*f;
        PP = PP + (P.Ts/N)*(dfdx*PP+PP*dfdx'+Q);
    end
    
    
    %If measurement has been received from sensor
    if(any(gps~=gps_d1))
        gps_d1=gps;
        
        dhdx = [...
                1,  0,  0,          0,              0,  0,  0;...
                0,  1,  0,          0,              0,  0,  0;...
                0,  0,  1,          0,              0,  0,  0;...
                0,  0,  0,          1,              0,  0,  0;...
                0,  0,  -cos(chi),  Vg*sin(chi),    1,  0,  -Va*sin(psi);...
                0,  0,  -sin(chi),  -Vg*cos(chi),   0,  1,  Va*cos(psi)];

        %wrap GPS course
        while (y_gps_course - chi)>pi, y_gps_course = y_gps_course - 2*pi; end
        while (y_gps_course - chi)>pi, y_gps_course = y_gps_course + 2*pi; end
        
        R = diag([P.sig_gps_n^2, P.sig_gps_e^2, P.sig_gps_vg^2, (P.sig_gps_vg/Vg)^2, .3^2, .3^2, .1^2]);
        y = [y_gps_n; y_gps_e; y_gps_Vg; y_gps_course; Va*cos(psi)+wn-Vg*cos(chi); Va*sin(psi)+we-Vg*sin(chi)];
        for i=1:6
            Ci = dhdx(i,:);
            Li = PP*Ci'/(R(i,i)+Ci*PP*Ci');
            PP = (eye(7)-Li*Ci)*PP;
            xhat = xhat+Li*(y(i)-xhat(i));
        end
    end
end
