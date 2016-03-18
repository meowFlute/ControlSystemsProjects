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
    
    cphi = cos(phi);
    sphi = sin(phi);
    cth = cos(theta);
    sth = sin(theta);
    cpsi = cos(psi);
    spsi = sin(psi);
    
    Rv_v1 = [cpsi,spsi,0;
             -spsi,cpsi,0;
             0,0,1];
    Rv1_v2 = [cth,0,-sth;
             0,1,0;
             sth,0,cth];
    Rv2_b = [1,0,0;
             0,cphi,sphi;
             0,-sphi,cphi];
    Rv_b = Rv2_b * Rv1_v2 * Rv_v1;
    
    % gust is in the vehicle frame
    gust = Rv_b' * [u_wg;v_wg;w_wg];
    % Vwb is in the body frame
    Vwb = Rv_b * [w_ns;w_es;w_ds] + [u_wg;v_wg;w_wg];
    % airspeed vector in the body frame
    Vab = [u;v;w]-Vwb;
    
    % compute wind data in NED
    w_n = w_ns+gust(1);
    w_e = w_es+gust(2);
    w_d = w_ds+gust(3);
    
    % compute air data
    Va = norm(Vab);
    alpha = atan2(Vab(3),Vab(1));
    if Va == 0
        beta = 0;
    else
        beta = asin(Vab(2)/Va);
    end
    
    % compute things here
    ms = P.mass;
    g = P.gravity;
    CL0 = P.C_L_0;
    CLa = P.C_L_alpha;
    CLq = P.C_L_q;
    CLde = P.C_L_delta_e;
    Cd0 = P.C_D_0;
    Cda = P.C_D_alpha;
    Cdp = P.C_D_p;
    Cdq = P.C_D_q;
    Cdde = P.C_D_delta_e;
    Cy0 = P.C_Y_0;
    Cyb = P.C_Y_beta;
    Cyp = P.C_Y_p;
    Cyr = P.C_Y_r;
    Cyda = P.C_Y_delta_a;
    Cydr = P.C_Y_delta_r;
    
    % compute Cd
    AR = (P.b)^2/P.S_wing;
    Cd = Cdp + (CL0+CLa*alpha)^2/(pi()*P.e*AR);
    % compute Cl
    sigma = (1+exp(-P.M*(alpha-P.alpha0))+exp(P.M*(alpha+P.alpha0))) / ...
                ((1+exp(-P.M*(alpha-P.alpha0)))*(1+exp(P.M*(alpha+P.alpha0))));
    CL = (1-sigma)*(CL0+CLa*alpha)+sigma*(2*sign(alpha)*(sin(alpha))^2*cos(alpha));
    
    % compute other aero stuff
    Cx = -Cd*cos(alpha)+CL*sin(alpha);
    Cxq = -Cdq*cos(alpha)+CLq*sin(alpha);
    Cxde = -Cdde*cos(alpha)+CLde*sin(alpha);
    Cz = -Cd*sin(alpha)-CL*cos(alpha);
    Czq = -Cdq*sin(alpha)-CLq*cos(alpha);
    Czde = -Cdde*sin(alpha)-CLde*cos(alpha);
    
    rho = P.rho;
    
    % compute external forces and torques on aircraft
     
    if Va == 0
        bVa2 = 0;
        cVa2 = 0;
    else
        bVa2 = P.b/Va/2;
        cVa2 = P.c/Va/2;
    end
    
    Force(1) =  -ms*g*sth + (rho*Va^2*P.S_wing/2)*(Cx+Cxq*cVa2*q+Cxde*delta_e) + (rho*P.S_prop*P.C_prop/2)*((P.k_motor*delta_t)^2-Va^2);
    Force(2) =  ms*g*cth*sphi + (rho*Va^2*P.S_wing/2)*(Cy0 + Cyb*beta + Cyp*bVa2*p + Cyr*bVa2*r + Cyda*delta_a + Cydr*delta_r);
    Force(3) =  ms*g*cth*cphi + (rho*Va^2*P.S_wing/2)*(Cz + Czq*cVa2*q + Czde*delta_e);
    
    Cl0 = P.C_ell_0;
    Clb = P.C_ell_beta;
    Clp = P.C_ell_p;
    Clr = P.C_ell_r;
    Clda = P.C_ell_delta_a;
    Cldr = P.C_ell_delta_r;
    Cm0 = P.C_m_0;
    Cma = P.C_m_alpha;
    Cmq = P.C_m_q;
    Cmde = P.C_m_delta_e;
    Cn0 = P.C_n_0;
    Cnb = P.C_n_beta;
    Cnp = P.C_n_p;
    Cnr = P.C_n_r;
    Cnda = P.C_n_delta_a;
    Cndr = P.C_n_delta_r;
    komega = P.k_Omega;
    ktp = P.k_T_P;
    
    Torque(1) = (rho*Va^2*P.S_wing/2) * P.b * (Cl0 + Clb*beta + Clp*bVa2*p + Clr*bVa2*r + Clda*delta_a + Cldr*delta_r) + (-ktp*(komega*delta_t)^2);
    Torque(2) = (rho*Va^2*P.S_wing/2) * P.c * (Cm0 + Cma*alpha + Cmq*cVa2*q + Cmde*delta_e);   
    Torque(3) = (rho*Va^2*P.S_wing/2) * P.b * (Cn0 + Cnb*beta + Cnp*bVa2*p + Cnr*bVa2*r + Cnda*delta_a + Cndr*delta_r);
   
    out = [Force'; Torque'; Va; alpha; beta; w_n; w_e; w_d];
end



