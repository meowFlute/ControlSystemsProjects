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
    
    mass = P.mass;
    g = P.g;
    rho = P.rho;
    
    
    
    % vehicle to body frame
    Rvb = [ ...
        cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta); ...
        sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);...
        cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)...
        ];
    
    % body frame to vehicle frame
    Rbv = Rvb';
    
    %% Calculate Wind Data
    
    % convert gust from body frame to vehicle frame
    wGustBody = [u_wg; v_wg; w_wg];
    wGustVehicle = Rbv*wGustBody;
    
    % compute wind data in NED
    w_n = wGustVehicle(1)+w_ns;
    w_e = wGustVehicle(2)+w_es;
    w_d = wGustVehicle(3)+w_ds;
    
    %% compute Va in body frame
    
    % convert wind into body frame
    Vwb = Rvb*[w_ns; w_es; w_ds] + wGustBody;
    
    % calculate Va in the body frame
    Vab = [u; v; w]-Vwb;
    
    % compute air data
    Va = norm(Vab);
    alpha = atan2(Vab(3),Vab(1));
    beta = 0;
    if (Va ~= 0)
        beta = asin(Vab(2)/Va);
    end
    
    %% compute external forces and torques on aircraft
    
    AR = (P.b)^2/P.S_wing;
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
    
    % Force due to gravity
    Fg = [-mass*g*sin(theta);...
        mass*g*cos(theta)*sin(phi);...
        mass*g*cos(theta)*cos(phi)...
        ];
    
    %Force due to aerodynamic effects
    Faero = 0;
    if (Va ~= 0)
    
        Faero = .5*rho*Va^2*P.S_wing*   [Cx+Cxq*P.c/(2*Va)*q+Cxde*delta_e;...
                                         P.C_Y_0+P.C_Y_beta*beta+P.C_Y_p*P.b/(2*Va)*p+P.C_Y_r*P.b/(2*Va)*r+P.C_Y_delta_a*delta_a+P.C_Y_delta_r*delta_r;...
                                         Cz+Czq*P.c/(2*Va)*q+Czde*delta_e];
    end
    
    
    %Force due to the propeller
    Fprop = .5*rho*P.S_prop*P.C_prop*[(P.k_motor*delta_t)^2-Va^2;0;0];
    
    %Sum of Forces;
    Force = Fg+Faero+Fprop;
    
    %Torque
    T1 = 0;
    if (Va ~= 0)
    T1 = .5*rho*Va^2*P.S_wing*[P.b*(P.C_ell_0+P.C_ell_beta*beta+P.C_ell_p*P.b/(2*Va)*p+P.C_ell_r*P.b/(2*Va)*r+P.C_ell_delta_a*delta_a+P.C_ell_delta_r*delta_r);...
                               P.c*(P.C_m_0+P.C_m_alpha*alpha+P.C_m_q*P.c/(2*Va)*q+P.C_m_delta_e*delta_e);...
                               P.b*(P.C_n_0+P.C_n_beta*beta+P.C_n_p*P.b/(2*Va)*p+P.C_n_r*P.b/(2*Va)*r+P.C_n_delta_a*delta_a+P.C_n_delta_r*delta_r)];
    end
    T2 = [-P.k_T_P*(P.k_Omega*delta_t)^2;0;0];
    
    Torque = T1+T2;

    out = [Force; Torque; Va; alpha; beta; w_n; w_e; w_d];
end