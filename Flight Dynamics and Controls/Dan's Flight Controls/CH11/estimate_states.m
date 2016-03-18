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
%   thetahat - estimated pitch angle, 
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
   
    persistent alpha alphaGPS;
    alpha = exp(-225*P.Ts);
    alphaGPS = exp(-1.5*P.Ts_gps);
    persistent phat qhat rhat pnhat pehat chihat Vghat;
    persistent phihat thetahat;
    persistent y_static_presd1 y_diff_presd1;
    persistent y_accel_xd1 y_accel_yd1 y_accel_zd1;
    persistent wnhat wehat psihat;
    
    persistent P1 Q1 R1 P2 Q2 R2;
    
    if t == 0
        % EKF 1 Values
        phat = y_gyro_x;
        qhat = y_gyro_y;
        rhat = y_gyro_z;
        hhat = y_static_presd1 / P.rho / P.gravity;
        Vahat = sqrt(2 / P.rho * y_diff_presd1);
        
        phihat = 0;%atan(y_accel_yd1 / y_accel_zd1);
        thetahat = 0;%asin(y_accel_xd1 / P.gravity);        

        y_static_presd1 = 0;
        y_diff_presd1 = 0;
        
        % EKF 2 Values
        pnhat = y_gps_n;
        pehat = y_gps_e;
        chihat = y_gps_course;
        Vghat = y_gps_Vg;
        
        wnhat = P.wind_n;
        wehat = P.wind_e;
        psihat = 0;
        
        %%% P,Q,R %%%
        P1 = diag([15*pi()/180,...   % variance in phi measurement
                   15*pi()/180]);    % variance in theta measurement
               
        Q1 = diag([0.0000001;       % process noise
                   0.0000001]);
        R1 = diag([P.sigma_accel^2;     % measurement of Ax
                   P.sigma_accel^2;     % measurement of Ay
                   P.sigma_accel^2;]*100);  % measurement of Az
               
        P2 = diag([3^2;     % variance in position north
                   3^2;     % variance in position east
                   0.1^2;   % variance in Vg
                   0.1^2;   % variance in chi
                   0.5^2;   % variance in wind north
                   0.5^2;   % variance in wind east
                   0.1^2]); % variance in psi
        Q2 = diag([0.1;0.1;0.1;0.1;0.1;0.1;0.1]);
        R2 = diag([P.sigma_n_gps^2;
                   P.sigma_e_gps^2;
                   1^2;
                   1^2;
                   1^2;
                   1^2]);
    end
    
    phat = LPF(y_gyro_x,phat,alpha);
    qhat = LPF(y_gyro_y,qhat,alpha);
    rhat = LPF(y_gyro_z,rhat,alpha);    
    y_static_presd1 = LPF(y_static_pres,y_static_presd1,alpha);
    y_diff_presd1 = LPF(y_diff_pres,y_diff_presd1,alpha);
    hhat = y_static_presd1 / P.rho / P.gravity;
    Vahat = sqrt(2 / P.rho * y_diff_presd1);
    
    % EKF #1: estimating [phi;theta]
    x0 = [phihat;thetahat];
    u = [y_gyro_x;y_gyro_y;y_gyro_z;Vahat];
    yn = [y_accel_x;y_accel_y;y_accel_z];
    STR = EKF1(x0,u,yn,P.Ts,10,P1,Q1,R1,P.gravity);
    P1 = STR.P;
    xhat = STR.x;
    phihat = xhat(1);
    thetahat = xhat(2);
    
    pnhat = LPF(y_gps_n,pnhat,alphaGPS);
    pehat = LPF(y_gps_e,pehat,alphaGPS);
    chihat = LPF(y_gps_course,chihat,alphaGPS);
    Vghat = LPF(y_gps_Vg,Vghat,alphaGPS);
    psihat = chihat;
    
    % EKF #2: estimating [pn;pe;Vg;chi;wind_n;wind_e;psi]
    x0 = [pnhat;pehat;Vghat;chihat;wnhat;wehat;psihat];
    u = [Vahat;qhat;rhat;phihat;thetahat];
    yn = [y_gps_n;y_gps_e;y_gps_Vg;y_gps_course];
    STR = EKF2(x0,u,yn,P.Ts,10,P2,Q2,R2,P.gravity);
    xhat = STR.x;
    P2 = STR.P;
    pnhat = xhat(1);
    pehat = xhat(2);
    Vghat = xhat(3);
    chihat = xhat(4);
    wnhat = xhat(5);
    wehat = xhat(6);
    psihat = xhat(7);
    
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
%{
function y = clamp(x)
    if x > 1, y = 1;
    elseif x < -1, y = -1;
    else y = x;
    end
end
%}
function yn1 = LPF(u,y,alpha)
    yn1 = (1.0-alpha) * y + alpha * u;
end

function S = EKF2(x0,u,y,Ts,N,P,Q,R,g)
    xhat = x0;
    ts = Ts / N;
    
    % input values 
    Va = u(1);
    q = u(2);
    r = u(3);
    phi = u(4);
    theta = u(5);
    
    % advance state using prediction
    for i=1:N
        pn = xhat(1);
        pe = xhat(2);
        Vg = xhat(3);
        chi = xhat(4);
        wn = xhat(5);
        we = xhat(6);
        psi = xhat(7);
        psidot = q*(sin(phi)/cos(theta))+r*(cos(phi)/cos(theta));
        
        f = [Vg*cos(theta);
             Vg*sin(theta);
             ((Va*cos(psi)+wn)*(-Va*psidot*sin(psi))+(Va*sin(psi)+we)*(Va*psidot*cos(psi))) / Vg;
             (g/Vg)*tan(phi)*cos(chi-psi);
             0;
             0;
             psidot];
        diff = ts*f;
        xhat = xhat + diff;
        pn = xhat(1);
        pe = xhat(2);
        Vg = xhat(3);
        chi = xhat(4);
        wn = xhat(5);
        we = xhat(6);
        psi = xhat(7);
        
        A = [0,0,cos(chi),-Vg*sin(chi),0,0,0;
             0,0,sin(chi),Vg*cos(chi),0,0,0;
             0,0,-((Va*cos(phi)+wn)*(-Va*psidot*sin(psi))+(Va*sin(psi)+we)*(Va*psidot*cos(psi)))/Vg^2,...
                0,-psidot*Va*sin(psi),psidot*Va*cos(psi),-psidot*Va*(wn*cos(psi)+we*sin(psi))/Vg;
             0,0,-g/Vg^2*tan(phi)*cos(chi-psi),-g/Vg*tan(phi)*sin(chi-psi),0,0,g/Vg*tan(phi)*sin(chi-psi);
             0,0,0,0,0,0,0;
             0,0,0,0,0,0,0;
             0,0,0,0,0,0,0];
        P = P * ts * (A*P+P*A'+Q);
    end
    
    I = eye(7);
    gps = [y(1);%gps_n;
           y(2);%gps_e;
           y(3);%gps_Vg;
           y(4);%gps_chi;
           0;
           0];
    for i=1:1
        pn = xhat(1);
        pe = xhat(2);
        Vg = xhat(3);
        chi = xhat(4);
        wn = xhat(5);
        we = xhat(6);
        psi = xhat(7);
        
        C = [1,0,0,0,0,0,0;
             0,1,0,0,0,0,0;
             0,0,1,0,0,0,0;
             0,0,0,1,0,0,0;
             0,0,-cos(chi),Vg*sin(chi),1,0,-Va*sin(psi);
             0,0,-sin(chi),-Vg*cos(chi),0,1,Va*cos(psi)];
        
        L = (P*C')*inv(C*P*C'+R);
        P = (I-L*C)*P;
        h = [pn;
             pe;
             Vg;
             chi;
             Va*cos(phi)+wn-Vg*cos(chi);
             Va*sin(phi)+we-Vg*sin(chi)];
        addition = L*(gps-h);
        xhat = xhat + addition;
    end
    
    S.x = xhat;
    S.P = P;
end

function S = EKF1(x0,u,y,Ts,N,P,Q,R,g)
    xhat = x0;
    p = u(1);
    q = u(2);
    r = u(3);
    Va = u(4);
    accel = y;
    ts = Ts / N;
    for i=1:N
        phi = xhat(1); % get new phi
        theta = xhat(2); % get new theta
        f = [p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);
             q*cos(phi)-r*sin(phi)];
        xhat = xhat + ts * f;
        
        phi = xhat(1); % get new phi
        theta = xhat(2); % get new theta
        A = [q*cos(phi)*tan(theta)-r*sin(phi)*tan(theta), (q*sin(phi)-r*cos(phi))/(cos(phi))^2;
             -q*sin(phi)-r*cos(phi), 0];
        P = P + ts * (A*P + P*A' + Q);
    end
    
    I = eye(2);
    for i=1:3
        phi = xhat(1); % get new phi
        theta = xhat(2); % get new theta
        
        C = [0,                      (q*Va+g)*cos(theta);
             -g*cos(phi)*cos(theta), -r*Va*sin(theta)-p*Va*cos(theta)+g*sin(phi)*sin(theta);
             g*sin(phi)*cos(theta),  (q*Va + g*cos(phi))*sin(theta)];
        Ci = C(i,:);
        Ri = R(i,i);
        L = (P*Ci') / (Ci*P*Ci'+Ri);
        P = (I-L*Ci)*P;
        
        phi = xhat(1); % get new phi
        theta = xhat(2); % get new theta
        h = [q*Va*sin(theta)+g*sin(theta);
             r*Va*cos(theta)-p*Va*sin(theta)-g*cos(theta)*sin(phi);
             -q*Va*cos(theta)-g*cos(theta)*cos(phi)];
        addition = L * (accel(i) - h(i));
        if(abs((accel(i) - h(i))) < 1)
            xhat = xhat + addition;
        end
    end
    
    S.x = xhat;
    S.P = P;
end
