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
   persistent gpsn1 gpse1 gpsh1 gpsv1 gpscourse1
   if(isempty(gyrox1))
      gyrox1 = y_gyro_x;
      gyroy1 = y_gyro_y;
      gyroz1 = y_gyro_z;
      staticpres1 = y_static_pres;
      diffpres1 = y_diff_pres;
      accelx1 = y_accel_x;
      accely1 = y_accel_y;
      accelz1 = y_accel_z;
      gpsn1 = y_gps_n;
      gpse1 = y_gps_e;
      gpsh1 = y_gps_h;
      gpsv1 = y_gps_Vg;
      gpscourse1 = y_gps_course;
   end
   phat = LPF(a, P.Ts, y_gyro_x, gyrox1);
   qhat = LPF(a, P.Ts, y_gyro_y, gyroy1);
   rhat = LPF(a, P.Ts, y_gyro_z, gyroz1);
   hhat = LPF(a, P.Ts, y_static_pres, staticpres1)/P.rho/P.g;
   Vahat = sqrt(2/P.rho*LPF(a, P.Ts, y_diff_pres, diffpres1));
   

   phihat = atan(LPF(a, P.Ts, y_accel_y, accely1)/LPF(a, P.Ts, y_accel_z, accelz1));
   thetahat = asin(LPF(a, P.Ts, y_accel_x, accelx1)/P.g);
   
    pnhat = LPF(agps, P.Ts_gps, gpsn1, y_gps_n);
    pehat = LPF(agps, P.Ts_gps, gpse1, y_gps_e);
    chihat = LPF(agps, P.Ts_gps, gpscourse1, y_gps_course);
    Vghat = LPF(agps, P.Ts_gps, gpsv1, y_gps_Vg);
    
    psihat = chihat;
    wnhat = 0;
    wehat = 0;
    %Kalman Filter #1
    
    dfdx1 = 
    %Kalman Filter #2
    psidothat = qhat*sin(phihat)/cos(thetahat)+rhat*cos(phihat)/cos(thetahat);
    vgdothat = (Vahat*cos(psihat)+wnhat)*(-Vahat*psidothat*sin(psihat))+(Vahat*sin(psihat)+wehat)*(Vahat*psidothat*cos(psihat))/Vghat;
    dfdx = [0   0   cos(chihat) -Vghat*sin(chihat)  0   0   0;...
            0   0   sin(chihat) Vghat*cos(chihat)   0   0   0;...
            0   0   -vgdothat/Vgdot     0   -psidothat*Vahat*sin(psihat)    psidothat*Vahat*cos(psihat) -psidothat*Vahat*(wnhat*cos(psihat)+
    
    wnhat = 0;
    wehat = 0;
    psihat = 0;
   
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

function x = EKF1(phi0, theta0, p, q, r, Va, Ax, Ay, Az, Tout, N, PP, Q)

    xhat = [phi0, theta0];
    u = [p, q, r, Va];
    for i = 0:N
       xhat = xhat + (Tout/N)*f(xhat, u);
       A = [q*cos(phi)*tan(theta)-r*sin(phi)*tan(theta) (q*sin(phi)-r*cos(phi))/cos(theta)^2;...
            -q*sin(phi)-r*cos(phi)  0];
       PP = PP + (Tout/N)*(A*PP+PP*A'+Q);
    end
    
end
function y = f1(x,u)
    phi = x(1);
    theta = x(2);
    p = u(1);
	q = u(2);
	r = u(3);
	Va = u(4);
    
    y = [p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);...
        q*cos(phi)-r*sin(phi)];
end
function y = h1(x,u, P)
    phi = x(1);
    theta = x(2);
    p = u(1);
	q = u(2);
	r = u(3);
	Va = u(4);
    
    y = [q*Va*sin(theta)+P.g*sin(theta);...
        r*Va*cos(theta)-p*Va*sin(theta)-P.g*cos(theta)*sin(phi);...
        -q*Va*cos(theta)-g*cos(theta)*cos(phi)];
end
